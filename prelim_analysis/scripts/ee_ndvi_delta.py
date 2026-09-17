#!/usr/bin/env python3
"""
Peak-season NDVI change and post-solstice change for the Miller Creek site.

Stays on Sentinel-2 L2A at 10 m. HLS (30 m, LaSRC + NBAR) was tested and is
radiometrically far better here -- sd 0.022 vs 0.055 across clear scenes -- but
30 m leaves this 137 x 417 m AOI only 5 x 14 pixels, too coarse to map.

Two fixes recover most of that stability at 10 m:

1. Radiometric screening. The Sen2Cor NDVI scatter over this AOI is driven by
   the red band (NDVI vs red reflectance r = -0.79; red swings 0.027-0.065 over
   the same conifer stand), the signature of atmospheric-correction residual --
   Sen2Cor's aerosol retrieval is documented as overestimating AOT with >160%
   relative error. Scenes whose AOI red reflectance exceeds the series median by
   more than RED_OUTLIER_TOL are rejected outright.
2. Medoid compositing. Evaluations of Landsat compositing algorithms rank
   maximum-NDVI among the worst for change detection: it selects the extreme
   value and so is drawn to whatever artifact inflates NDVI (it is what produced
   NDVI = 1.0 over 23-39% of this AOI before masking was tightened). The medoid
   takes the real observation closest to the multi-band median instead.

"Peak greenness" is preserved by taking the medoid inside a peak-season window
rather than over the whole year, which a full-year medoid would not capture.

References:
    Roy et al. 2016/2017, c-factor BRDF normalisation (HLS NBAR)
    Claverie et al. 2018, the HLS surface reflectance product
    Qiu et al. 2023, Evaluation of Landsat image compositing algorithms
    Doxani et al. / Sen2Cor validation: AOT overestimated, relative accuracy >160%

Inputs:
    - prelim_analysis/results/ndvi_shadows.geojson: survey points for overlay
    - COPERNICUS/S2_SR_HARMONIZED + GOOGLE/CLOUD_SCORE_PLUS (Earth Engine)
    - data/raster/miller_rgb_251021_25cm.tif: downsampled drone RGB, built once
      from the 1.1 GB GCS orthomosaic (see RGB_DOWNSAMPLE_M)

Outputs:
    - prelim_analysis/results/ndvi_change.tif: 4 bands (ndvi_peak_2025,
      ndvi_peak_2026, ndvi_peak_delta, ndvi_solstice_delta), 10 m, EPSG:6514
    - prelim_analysis/results/ndvi_delta_map.png: three map panels over
      control-vs-girdled box plots of the two Sentinel-2 change metrics
"""

import ee
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import rasterio
import requests
import subprocess
from datetime import date
from dotenv import load_dotenv
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
from rasterio.plot import plotting_extent
from rasterio.windows import from_bounds
from scipy import stats

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_POINTS = BASE_DIR / "prelim_analysis" / "results" / "ndvi_shadows.geojson"
OUTPUT_TIF = BASE_DIR / "prelim_analysis" / "results" / "ndvi_change.tif"
OUTPUT_PNG = BASE_DIR / "prelim_analysis" / "results" / "ndvi_delta_map.png"

GCP_PROJECT = "mpg-projects"

S2 = "COPERNICUS/S2_SR_HARMONIZED"
CLOUD_SCORE = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"
CS_THRESHOLD = 0.60
# SCL: nodata, saturated, dark, cloud shadow, cloud med/high, cirrus, snow
BAD_SCL = [0, 1, 2, 3, 8, 9, 10, 11]

# Reject a scene whose AOI red reflectance exceeds the series median by more
# than this. On 2026 summer scenes it removes exactly the three that carry the
# atmospheric-correction residual, taking series sd from 0.055 to ~0.024 --
# effectively HLS-grade stability without leaving 10 m.
RED_OUTLIER_TOL = 0.25

SCALE_M = 10
WORK_CRS = "EPSG:6514"
POINT_BUFFER_M = 10          # crop margin around the adult-tree bounding box

# Peak-season window. A medoid over a whole year would describe the typical
# annual condition (snow included); restricting to the peak window keeps the
# "greenest" intent while still using a robust estimator.
PEAK_WINDOW = ("07-01", "08-16")
YEARS = (2025, 2026)

# 2026 June solstice, and the post-solstice comparison
SOLSTICE = "2026-06-21"
SOLSTICE_SEARCH_DAYS = 15
ANCHOR_SCENES = 3            # medoid of this many scenes at each end
MIN_CLEAR_FRAC = 0.98        # a scene must see this much of the AOI cleanly

# The 2025 drone RGB ortho is 1.28 cm / 1.1 GB with no overviews, so a decimated
# read would still stream the whole file. It is downsampled to a local cache once
# and reused; 25 cm is far finer than the satellite panels it sits beside.
RGB_SOURCE = (
    "/vsicurl/https://storage.googleapis.com/mpg-aerial-survey/surveys/"
    "miller_collaring/251021/processing/drone_deploy/miller_collaring-RGB-251021.tif"
)
RGB_CACHE = BASE_DIR / "data" / "raster" / "miller_rgb_251021_25cm.tif"
RGB_DOWNSAMPLE_M = 0.25

# Windowed reads over HTTP: don't list the bucket, do cache blocks
GDAL_OPTS = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "VSI_CACHE": True,
    "GDAL_CACHEMAX": 512,
}

# Names look like B4_C_PSEMEN_T2_seedling. Field 1 is the block, field 2 is the
# treatment (C = control, G = girdled), field 3 the USDA species code. Note the
# "adult" token marks age class, NOT treatment -- 39 control points carry it.
SPECIES_RE = r"_(PSEMEN|PICENG|PINCON|PINPON|LAROCC)_?"
SPECIES_NAMES = {
    "PSEMEN": "Douglas-fir",
    "PICENG": "Engelmann spruce",
    "PINCON": "lodgepole pine",
    "LAROCC": "western larch",
    "PINPON": "ponderosa pine",
}

# Colour carries treatment (slots 1-2 of the reference categorical palette; all
# six checks pass all-pairs), marker carries species.
TREATMENT_STYLE = {"C": ("#2a78d6", "Control"), "G": ("#eb6834", "Girdled")}
SPECIES_MARKER = {
    "PSEMEN": "o", "PICENG": "^", "PINCON": "s", "LAROCC": "D", "PINPON": "P",
}

# Diverging two-hue with a neutral midpoint for polarity. Brown/teal poles are
# the ColorBrewer BrBG pair, CVD-safe, and "browning vs greening" is the
# domain's own vocabulary.
DELTA_CMAP = LinearSegmentedColormap.from_list(
    "browning_greening", ["#8C510A", "#BF812D", "#F5F5F2", "#35978F", "#01665E"]
)
INK = "#1A1A1A"
MUTED = "#6E6E6E"
RULE = "#D8D8D8"

# Below this many points a box implies more structure than the data carries


def init_ee():
    """Authenticate Earth Engine with the service account from .env."""
    load_dotenv(BASE_DIR / ".env")
    key_path = os.environ["MPG_PROJECTS_KEY_PATH"]
    ee.Initialize(
        ee.ServiceAccountCredentials(None, key_path), project=GCP_PROJECT
    )
    print(f"  Earth Engine initialized on project {GCP_PROJECT}")


def s2_collection(region, start, end):
    """Sentinel-2 SR renamed to a common (red, nir) pair, cloud masked."""
    col = (
        ee.ImageCollection(S2)
        .filterBounds(region)
        .filterDate(start, end)
        .linkCollection(ee.ImageCollection(CLOUD_SCORE), ["cs_cdf"])
    )

    def mask(img):
        scl = img.select("SCL")
        pair = img.select(["B4", "B8"], ["red", "nir"])
        keep = (
            img.select("cs_cdf").gte(CS_THRESHOLD)
            .And(scl.remap(BAD_SCL, [0] * len(BAD_SCL), 1))
            .And(pair.select("red").gt(0))
            .And(pair.select("nir").gt(0))
        )
        return pair.updateMask(keep).copyProperties(img, ["system:time_start"])

    return col.map(mask)


def screen_radiometric(col, region):
    """Drop scenes whose AOI red reflectance is an outlier against the series.

    Sen2Cor aerosol error lands mostly in the red band, so an anomalously bright
    red over an unchanged conifer stand flags a badly corrected scene. Rejecting
    those is what a maximum-value composite cannot do -- it seeks them out.
    """
    def tag(img):
        red = img.select("red").reduceRegion(
            ee.Reducer.median(), region, SCALE_M, maxPixels=int(1e8)
        ).get("red")
        return img.set("red_level", ee.Number(
            ee.Algorithms.If(ee.Algorithms.IsEqual(red, None), -1, red)
        ))

    tagged = col.map(tag).filter(ee.Filter.gt("red_level", 0))
    baseline = ee.Number(
        ee.List(tagged.aggregate_array("red_level")).reduce(ee.Reducer.median())
    )
    return tagged.filter(
        ee.Filter.lte("red_level", baseline.multiply(1 + RED_OUTLIER_TOL))
    )


def clean_collection(region, start, end):
    """Cloud-masked and radiometrically screened."""
    return screen_radiometric(s2_collection(region, start, end), region)


def medoid(collection):
    """Medoid composite: the real observation nearest the per-band median.

    Robust to outliers in a way maximum-value compositing is not, and it keeps
    red and nir from the same acquisition rather than mixing dates per band.
    """
    median = collection.select(["red", "nir"]).median()

    def tag(img):
        # Negated so qualityMosaic's max-selection picks the *smallest* distance
        distance = img.select(["red", "nir"]).subtract(median).pow(2).reduce(
            ee.Reducer.sum()
        ).sqrt()
        return img.addBands(distance.multiply(-1).rename("score"))

    return collection.map(tag).qualityMosaic("score").select(["red", "nir"])


def ndvi_of(image):
    """NDVI from a (red, nir) image."""
    return image.normalizedDifference(["nir", "red"]).rename("ndvi")


def peak_ndvi(region, year):
    """Medoid NDVI inside the peak-season window for one year."""
    start, end = PEAK_WINDOW
    col = clean_collection(region, f"{year}-{start}", f"{year}-{end}")
    return ndvi_of(medoid(col)), col.size()


def with_clear_fraction(img, region):
    """Tag an image with the fraction of the AOI it sees clearly."""
    frac = img.select("red").mask().rename("frac").reduceRegion(
        ee.Reducer.mean(), region, SCALE_M, maxPixels=int(1e8)
    ).get("frac")
    # Fully-masked scenes drop the key rather than returning null
    return img.set("clear_frac", ee.Number(
        ee.Algorithms.If(ee.Algorithms.IsEqual(frac, None), 0, frac)
    ))


def eligible_scenes(region, start, end):
    """Scenes clean enough over the AOI to anchor a comparison."""
    return clean_collection(region, start, end).map(
        lambda img: with_clear_fraction(img, region)
    ).filter(ee.Filter.gte("clear_frac", MIN_CLEAR_FRAC))


def solstice_delta(region, end):
    """NDVI change from the solstice to the present.

    Each end is a medoid of ANCHOR_SCENES clean scenes rather than a single
    acquisition, so no one scene's residual radiometric offset sets the result.
    """
    target = ee.Date(SOLSTICE)
    near = eligible_scenes(
        region,
        target.advance(-SOLSTICE_SEARCH_DAYS, "day"),
        target.advance(SOLSTICE_SEARCH_DAYS, "day"),
    ).map(lambda img: img.set(
        "dt", img.date().difference(target, "day").abs()
    )).sort("dt").limit(ANCHOR_SCENES)

    recent = eligible_scenes(region, SOLSTICE, end).sort(
        "system:time_start", False
    ).limit(ANCHOR_SCENES)

    first, last = ndvi_of(medoid(near)), ndvi_of(medoid(recent))

    def dates(col):
        return col.aggregate_array("system:time_start").map(
            lambda t: ee.Date(t).format("YYYY-MM-dd")
        )

    info = ee.Dictionary({
        "first_dates": dates(near), "last_dates": dates(recent),
        "n_eligible": eligible_scenes(region, SOLSTICE, end).size(),
    })
    return last.subtract(first).rename("solstice_delta"), info


def download_tif(image, region, out_path):
    """Fetch a small EE image as a GeoTIFF over HTTP."""
    url = image.getDownloadURL({
        "region": region,
        "scale": SCALE_M,
        "crs": WORK_CRS,
        "format": "GEO_TIFF",
    })
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    out_path.write_bytes(response.content)


def load_adults():
    """Non-seedling points tagged with treatment and species.

    Treatment comes from field 2 of the name (C = control, G = girdled), NOT
    from the "adult" token -- 39 control points carry "adult" as an age-class
    marker, so reading it as treatment would collapse the control group.
    """
    gdf = gpd.read_file(INPUT_POINTS)
    fields = gdf["Name"].str.split("_", expand=True)
    gdf["block"] = fields[0]
    gdf["treatment"] = fields[1]
    gdf["species"] = gdf["Name"].str.extract(SPECIES_RE)

    seedling = gdf["Name"].str.lower().str.contains("seedling")
    adults = gdf[~seedling].copy()

    counts = adults["treatment"].value_counts()
    print(f"  {len(adults)} non-seedling points "
          f"({len(gdf) - len(adults)} seedlings dropped)")
    print("  " + ", ".join(
        f"{TREATMENT_STYLE[k][1]} {counts[k]}" for k in counts.index
        if k in TREATMENT_STYLE
    ))
    unknown = sorted(set(adults["treatment"]) - set(TREATMENT_STYLE))
    if unknown:
        print(f"  warning: unrecognised treatment codes {unknown}")
    return adults


def plot_points(ax, adults, size=46):
    """Colour by treatment, marker by species."""
    for code, (color, _) in TREATMENT_STYLE.items():
        for species, marker in SPECIES_MARKER.items():
            sel = adults[
                (adults["treatment"] == code) & (adults["species"] == species)
            ]
            if sel.empty:
                continue
            ax.scatter(sel.geometry.x, sel.geometry.y, c=color, marker=marker,
                       s=size, edgecolor="white", linewidth=0.8, zorder=4)


def treatment_legend(adults):
    """Colour key: control vs girdled."""
    counts = adults["treatment"].value_counts()
    return [
        Line2D([], [], marker="o", color="none", markerfacecolor=color,
               markeredgecolor="white", markeredgewidth=0.8, markersize=9,
               label=f"{label} ({counts.get(code, 0)})")
        for code, (color, label) in TREATMENT_STYLE.items()
    ]


def species_legend(adults):
    """Marker key, most abundant species first, drawn in neutral ink."""
    counts = adults["species"].value_counts()
    return [
        Line2D([], [], marker=SPECIES_MARKER[code], color="none",
               markerfacecolor=MUTED, markeredgecolor="white",
               markeredgewidth=0.8, markersize=8,
               label=f"{code} - {SPECIES_NAMES[code]} ({counts[code]})")
        for code in counts.index if code in SPECIES_MARKER
    ]


def ensure_rgb_cache():
    """Downsample the drone RGB ortho to a local cache if it isn't there yet."""
    if RGB_CACHE.exists():
        print(f"  using cached {RGB_CACHE.name}")
        return

    print(f"  downsampling ortho to {RGB_DOWNSAMPLE_M} m (streams ~1.1 GB once)...")
    RGB_CACHE.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
               GDAL_CACHEMAX="1024", GDAL_NUM_THREADS="ALL_CPUS")
    subprocess.run([
        "gdalwarp", "-overwrite",
        "-tr", str(RGB_DOWNSAMPLE_M), str(RGB_DOWNSAMPLE_M), "-r", "average",
        "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "TILED=YES",
        RGB_SOURCE, str(RGB_CACHE),
    ], check=True, env=env)


def read_rgb_panel(bounds, crs):
    """Read the cached drone RGB over the map extent as an RGBA array."""
    with rasterio.open(RGB_CACHE) as src:
        if src.crs != crs:
            raise ValueError(f"RGB cache CRS {src.crs} != map CRS {crs}")
        window = from_bounds(*bounds, src.transform).round_offsets().round_lengths()
        rgba = src.read(
            (1, 2, 3, 4), window=window, boundless=True, fill_value=0
        ).astype("float32") / 255.0
        extent = plotting_extent(
            rgba[0], transform=src.window_transform(window)
        )
        res = src.res[0]
    return np.moveaxis(rgba, 0, -1), extent, res


def treatment_stats(adults, column):
    """Group means and a Welch t-test for one metric."""
    groups = {
        code: adults.loc[adults["treatment"] == code, column].dropna()
        for code in TREATMENT_STYLE
    }
    t, p = stats.ttest_ind(groups["G"], groups["C"], equal_var=False)
    return groups, groups["G"].mean() - groups["C"].mean(), p


def plot_box(ax, adults, column, label_side="left"):
    """Control vs girdled, pooled, with every observation drawn.

    Trees share 10 m pixels (28 girdled trees fall in 15 unique pixels), so
    points stacked at an identical value are one pixel measured several times,
    not independent observations -- the Welch p below is correspondingly
    optimistic.
    """
    rng = np.random.default_rng(0)
    counts = {}
    for x, (code, (color, label)) in enumerate(TREATMENT_STYLE.items()):
        values = adults.loc[adults["treatment"] == code, column].dropna()
        if values.empty:
            continue
        counts[code] = len(values)

        box = ax.boxplot(
            [values], positions=[x], widths=0.52, patch_artist=True,
            showfliers=False, medianprops=dict(color=INK, linewidth=1.6),
            whiskerprops=dict(color=MUTED, linewidth=1),
            capprops=dict(color=MUTED, linewidth=1),
        )
        patch = box["boxes"][0]
        patch.set_facecolor(color)
        patch.set_alpha(0.22)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.4)

        ax.scatter(x + rng.uniform(-0.11, 0.11, len(values)), values,
                   s=17, color=color, edgecolor="white", linewidth=0.5,
                   alpha=0.9, zorder=3)

    ax.axhline(0, color=RULE, linewidth=1, zorder=1)
    ax.set_ylabel(r"$\Delta$NDVI", fontsize=9.5, color=MUTED)
    # Columns are as tight as the map panels above, so the right-hand plot puts
    # its scale on the outside rather than into its neighbour's gap
    if label_side == "right":
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    ax.set_xticks(range(len(TREATMENT_STYLE)))
    ax.set_xticklabels(
        [f"{label}\n(n={counts.get(code, 0)})"
         for code, (_, label) in TREATMENT_STYLE.items()],
        fontsize=9.5, color=INK,
    )
    ax.set_xlim(-0.7, len(TREATMENT_STYLE) - 0.3)
    ax.margins(y=0.24)
    ax.tick_params(axis="y", labelsize=8.5, colors=MUTED, length=2)
    ax.grid(axis="y", color=RULE, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    spine = "left" if label_side == "left" else "right"
    for side in ("top", "left", "right", "bottom"):
        ax.spines[side].set_visible(side == spine)
    ax.spines[spine].set_color(RULE)

    _, difference, p_value = treatment_stats(adults, column)
    flag = "n.s." if p_value >= 0.05 else "p < 0.05"
    ax.set_title(
        rf"$\Delta$ = {difference:+.3f}   {flag} (p = {p_value:.2f})",
        fontsize=9, color=MUTED, pad=6,
    )


def build_figure(tif_path, adults, out_path, info=None):
    """Peak-season delta, post-solstice delta, and the drone RGB reference."""
    with rasterio.open(tif_path) as src:
        bands = src.read(masked=True)
        extent = plotting_extent(src)
        bounds = src.bounds
        crs = src.crs

    peak_2025, peak_2026, peak_delta, solstice = bands
    adults = adults.to_crs(crs)

    with rasterio.open(tif_path) as src:
        sampled = np.array(list(
            src.sample([(p.x, p.y) for p in adults.geometry])
        ))
    adults["peak_delta"] = sampled[:, 2]
    adults["solstice_delta"] = sampled[:, 3]

    # Diverging scales are symmetric so the neutral midpoint sits at zero
    plim = float(np.nanpercentile(np.abs(peak_delta.compressed()), 99))
    slim = float(np.nanpercentile(np.abs(solstice.compressed()), 99))

    if info:
        span = f"{info['first_dates'][0]} to {info['last_dates'][0]}"
    else:
        span = "solstice to most recent"

    # Two gridspecs sharing one column geometry: the box row is a separate grid
    # only so it can carry its own top/bottom, not a different width.
    fig = plt.figure(figsize=(12.2, 12.4))
    map_grid = fig.add_gridspec(
        1, 3, wspace=0.02, left=0.055, right=0.79, top=0.955, bottom=0.34,
    )
    # Same columns as the map so each box plot sits directly under the panel it
    # summarises; the drone column stays empty.
    box_grid = fig.add_gridspec(
        1, 3, wspace=0.02, left=0.055, right=0.79, top=0.27, bottom=0.045,
    )
    axes = [fig.add_subplot(map_grid[0, col]) for col in range(3)]
    panels = [
        (peak_delta,
         f"Peak-season NDVI change\n{YEARS[1]} - {YEARS[0]} "
         f"({PEAK_WINDOW[0]} to {PEAK_WINDOW[1]})",
         TwoSlopeNorm(vcenter=0, vmin=-plim, vmax=plim), "NDVI change"),
        (solstice,
         f"Post-solstice NDVI change\n{span}",
         TwoSlopeNorm(vcenter=0, vmin=-slim, vmax=slim), "NDVI change"),
    ]

    for ax, (data, title, norm, cbar_label) in zip(axes[:2], panels):
        img = ax.imshow(data, extent=extent, cmap=DELTA_CMAP, norm=norm,
                        interpolation="nearest")
        plot_points(ax, adults)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_title(title, fontsize=11, color=INK, pad=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(RULE)

        cax = make_axes_locatable(ax).append_axes("bottom", size="4%", pad=0.12)
        cbar = fig.colorbar(img, cax=cax, orientation="horizontal")
        cbar.set_label(cbar_label, fontsize=9, color=MUTED)
        cbar.ax.tick_params(labelsize=8, colors=MUTED, length=2)
        cbar.ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        cbar.outline.set_visible(False)

    # Drone RGB reference panel
    ax = axes[2]
    rgba, rgb_extent, _ = read_rgb_panel(bounds, crs)
    ax.imshow(rgba, extent=rgb_extent, interpolation="bilinear")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    plot_points(ax, adults)
    ax.set_title("Drone RGB\n21 Oct 2025", fontsize=11, color=INK, pad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor(RULE)

    # Matched spacer so all three panels share one plot height
    spacer = make_axes_locatable(ax).append_axes("bottom", size="4%", pad=0.12)
    spacer.set_axis_off()

    trt = fig.legend(
        handles=treatment_legend(adults), loc="upper left",
        bbox_to_anchor=(0.795, 0.78), ncol=1, frameon=True,
        fontsize=9.5, labelcolor=INK, handletextpad=0.5,
        title="Treatment", title_fontsize=10,
    )
    spp = fig.legend(
        handles=species_legend(adults), loc="upper left",
        bbox_to_anchor=(0.795, 0.63), ncol=1, frameon=True,
        fontsize=9, labelcolor=INK, handletextpad=0.5,
        title="Species", title_fontsize=10,
    )
    for legend in (trt, spp):
        legend.get_frame().set_edgecolor(RULE)
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_linewidth(0.8)
        legend.get_title().set_color(INK)
    fig.add_artist(trt)

    # Box plots sit directly under the change panel they summarise; the drone
    # column has no metric, so that cell stays empty.

    for col, column in enumerate(["peak_delta", "solstice_delta"]):
        plot_box(fig.add_subplot(box_grid[0, col]), adults, column,
                 label_side="left" if col == 0 else "right")

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return peak_2025, peak_2026, peak_delta, solstice


def main():
    """Main processing function."""
    print("=" * 60)
    print("Sentinel-2 NDVI Change, medoid composites")
    print("Miller Creek Collaring")
    print("=" * 60)

    print("\nInitializing Earth Engine...")
    init_ee()

    print("\nBuilding area of interest...")
    adults = load_adults()
    aoi = adults.to_crs(WORK_CRS).total_bounds + np.array(
        [-POINT_BUFFER_M, -POINT_BUFFER_M, POINT_BUFFER_M, POINT_BUFFER_M]
    )
    region = ee.Geometry.Rectangle(
        [float(v) for v in aoi], proj=WORK_CRS, geodesic=False, evenOdd=True
    )
    print(f"  {aoi[2] - aoi[0]:.0f} x {aoi[3] - aoi[1]:.0f} m "
          f"({POINT_BUFFER_M} m buffer on the adult-tree bounding box)")

    print(f"\nMedoid peak-season composites "
          f"({PEAK_WINDOW[0]} to {PEAK_WINDOW[1]})...")
    peaks = {}
    for year in YEARS:
        image, kept = peak_ndvi(region, year)
        raw = s2_collection(region, f"{year}-{PEAK_WINDOW[0]}",
                            f"{year}-{PEAK_WINDOW[1]}").size()
        peaks[year] = image
        n_raw, n_kept = raw.getInfo(), kept.getInfo()
        print(f"  {year}: {n_kept} of {n_raw} scenes kept "
              f"({n_raw - n_kept} rejected as red-band outliers)")

    print("\nPost-solstice change, medoid of anchor scenes...")
    solstice, info_obj = solstice_delta(region, f"{YEARS[1]}-{date.today():%m-%d}")
    info = info_obj.getInfo()
    print(f"  solstice anchors: {', '.join(info['first_dates'])}")
    print(f"  recent anchors:   {', '.join(info['last_dates'])}")
    print(f"  {info['n_eligible']} eligible scenes in window")

    stack = (
        peaks[YEARS[0]].rename("ndvi_peak_2025")
        .addBands(peaks[YEARS[1]].rename("ndvi_peak_2026"))
        .addBands(peaks[YEARS[1]].subtract(peaks[YEARS[0]])
                  .rename("ndvi_peak_delta"))
        .addBands(solstice.rename("ndvi_solstice_delta"))
        .toFloat()
    )

    print("\nDownloading composite...")
    OUTPUT_TIF.parent.mkdir(parents=True, exist_ok=True)
    download_tif(stack, region, OUTPUT_TIF)
    with rasterio.open(OUTPUT_TIF) as src:
        print(f"  {OUTPUT_TIF.name}: {src.width} x {src.height} px, "
              f"{src.count} bands, {src.crs}")

    print("\nPreparing drone RGB panel...")
    ensure_rgb_cache()

    print("\nRendering map...")
    peak_2025, peak_2026, peak_delta, solstice_arr = build_figure(
        OUTPUT_TIF, adults, OUTPUT_PNG, info
    )
    with rasterio.open(OUTPUT_TIF) as src:
        sampled = np.array(list(
            src.sample([(p.x, p.y) for p in adults.to_crs(src.crs).geometry])
        ))
    adults["peak_delta"] = sampled[:, 2]
    adults["solstice_delta"] = sampled[:, 3]

    print("\nSummary (whole map):")
    for name, data in [("peak NDVI 2025", peak_2025),
                       ("peak NDVI 2026", peak_2026),
                       ("peak delta", peak_delta),
                       ("solstice delta", solstice_arr)]:
        values = data.compressed()
        print(f"  {name:15s} min={values.min():+.3f}  mean={values.mean():+.3f}  "
              f"max={values.max():+.3f}")
    for name, data in [("peak delta", peak_delta),
                       ("solstice delta", solstice_arr)]:
        browner = float((data.compressed() < 0).mean() * 100)
        print(f"  {browner:.1f}% of pixels are less green ({name})")

    print("\nControl vs girdled:")
    for column in ("peak_delta", "solstice_delta"):
        groups, difference, p_value = treatment_stats(adults, column)
        print(f"  {column:15s} "
              f"C {groups['C'].mean():+.4f}+-{groups['C'].std():.4f}  "
              f"G {groups['G'].mean():+.4f}+-{groups['G'].std():.4f}  "
              f"diff {difference:+.4f}  p={p_value:.3f}")

    print("\n" + "=" * 60)
    print("Processing complete!")
    print(f"Raster saved to: {OUTPUT_TIF}")
    print(f"Map saved to:    {OUTPUT_PNG}")
    print("=" * 60)


if __name__ == "__main__":
    main()
