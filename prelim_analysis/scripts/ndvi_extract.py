#!/usr/bin/env python3
"""
Extract NDVI, RGB thumbnails, and a shadow index for Miller Creek seedling points.

Reads the 2025-10-21 DroneDeploy orthomosaics straight from GCS over /vsicurl
(windowed reads only -- the full rasters are never downloaded).

Multispectral band order (per DroneDeploy export): R, G, NIR, RE, unused, alpha
Visible band order: R, G, B, alpha

RTK quality is carried through from the Emlid export so a weak fix can be told
apart from a genuinely low-NDVI point: gps_rms_m is the Lateral (horizontal) RMS
and tilt_deg is the pole tilt angle. Either is rendered !!!NONE!!! when absent.

Shadow index follows Rikimaru et al. (2002) Forest Canopy Density mapping:
    SI = [(1 - R)(1 - G)(1 - B)] ^ (1/3)
with each visible band scaled to 0-1. Higher values mean darker / more shadowed.

Inputs:
    - data/tabular/miller_seedlings.csv: Emlid RTK survey points (EPSG:4326)
    - <GCS>/miller_collaring-MS-251021.tif: 5-band + alpha MS orthomosaic
    - <GCS>/miller_collaring-RGB-251021.tif: RGB + alpha orthomosaic

Outputs:
    - data/raster/seedling_thumbnails/2025/{seedling_id}.png: visible chip over
      the 3m buffer, annotated with a 0.5m-diameter ring and the extracted values
    - prelim_analysis/results/ndvi_shadows.geojson: points with ndvi_0o25,
      ndvi_0o1, shadow_idx, and the RTK quality fields gps_rms_m / tilt_deg
"""

import geopandas as gpd
import numpy as np
import os
import pandas as pd
import rasterio
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from rasterio.windows import Window, from_bounds

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_POINTS = BASE_DIR / "data" / "tabular" / "miller_seedlings.csv"
OUTPUT_THUMBS = BASE_DIR / "data" / "raster" / "seedling_thumbnails" / "2025"
OUTPUT_GEOJSON = BASE_DIR / "prelim_analysis" / "results" / "ndvi_shadows.geojson"

GCS_PREFIX = (
    "/vsicurl/https://storage.googleapis.com/mpg-aerial-survey/surveys/"
    "miller_collaring/251021/processing/drone_deploy"
)
INPUT_MS = f"{GCS_PREFIX}/miller_collaring-MS-251021.tif"
INPUT_RGB = f"{GCS_PREFIX}/miller_collaring-RGB-251021.tif"

# Multispectral band indices (1-based, rasterio convention)
BAND_RED = 1
BAND_NIR = 3
BAND_MS_ALPHA = 6
BAND_RGB_ALPHA = 4

# Extraction radii in meters
NDVI_RADII = {"ndvi_0o25": 0.25, "ndvi_0o1": 0.1}
THUMB_RADIUS = 3.0

# Thumbnail annotation. The chip is upsampled by SUPERSAMPLE before anything is
# drawn on it, so the ring and the label text carry that many times the pixel
# detail of the native 1.28 cm imagery and stay crisp when the report scales them.
SUPERSAMPLE = 2
RING_DIAMETER = 0.5          # meters, open circle drawn at the survey point
RING_COLOR = (255, 0, 0, 110)
RING_WIDTH = 2
LABEL_BG = (0, 0, 0, 185)
LABEL_FG = (255, 255, 255, 255)
LABEL_PAD = 8
LABEL_FONT_SIZE = 15
LABEL_LINE_SPACING = 3
MISSING_LABEL = "!!!NONE!!!"

# Thumbnails are rendered in a process pool; GDAL datasets are per-worker
N_WORKERS = max(1, (os.cpu_count() or 2) - 1)

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
]

# Windowed reads over HTTP: don't list the bucket, do cache blocks
GDAL_OPTS = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "VSI_CACHE": True,
    "GDAL_CACHEMAX": 512,
}


def load_points():
    """Load Emlid survey points as a GeoDataFrame in EPSG:4326."""
    print("Loading seedling points...")

    df = pd.read_csv(INPUT_POINTS)
    df = df.dropna(subset=["Longitude", "Latitude"]).reset_index(drop=True)
    df["idx"] = df.index

    # RTK quality: Lateral RMS is the horizontal accuracy that decides whether
    # the small buffers actually land on the seedling
    df["gps_rms_m"] = df["Lateral RMS"]
    df["tilt_deg"] = df["Tilt angle"]

    n_no_tilt = int(df["tilt_deg"].isna().sum())
    if n_no_tilt:
        print(f"  {n_no_tilt} point(s) have no tilt angle -> {MISSING_LABEL}")

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["Longitude"], df["Latitude"]),
        crs="EPSG:4326",
    )

    print(f"  {len(gdf)} points loaded")
    return gdf


def fmt(value):
    """Format an extracted index value for the thumbnail label."""
    return MISSING_LABEL if value is None else f"{value:+.3f}"


def fmt_quality(value, unit, decimals=3):
    """Format an RTK quality value, flagging absent entries loudly."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return MISSING_LABEL
    return f"{value:.{decimals}f} {unit}"


def load_font(size=15):
    """First available TrueType font, falling back to PIL's bitmap default."""
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def circle_mask(src, x, y, radius):
    """Return (window, mask) selecting pixels whose centers fall within radius.

    Falls back to the single pixel containing the point when the radius is
    smaller than a pixel. Returns (None, None) if the circle is off the raster.
    """
    window = from_bounds(
        x - radius, y - radius, x + radius, y + radius, src.transform
    ).round_offsets().round_lengths()

    # Pad by a pixel so pixel centers near the edge of the circle aren't clipped
    window = Window(
        window.col_off - 1, window.row_off - 1, window.width + 2, window.height + 2
    )
    window = window.intersection(Window(0, 0, src.width, src.height))
    if window.width <= 0 or window.height <= 0:
        return None, None

    rows, cols = np.mgrid[
        window.row_off:window.row_off + window.height,
        window.col_off:window.col_off + window.width,
    ]
    # Pixel centers in map coordinates (transform.xy flattens, so reshape back)
    px, py = rasterio.transform.xy(src.transform, rows, cols)
    px = np.asarray(px).reshape(rows.shape)
    py = np.asarray(py).reshape(rows.shape)
    mask = (px - x) ** 2 + (py - y) ** 2 <= radius ** 2

    if not mask.any():
        # Radius is sub-pixel: take the pixel the point lands in
        row, col = src.index(x, y)
        if not (0 <= row < src.height and 0 <= col < src.width):
            return None, None
        return Window(col, row, 1, 1), np.ones((1, 1), dtype=bool)

    return window, mask


def extract_ndvi(src, x, y, radius):
    """Mean NDVI over a circular buffer, or None where no valid pixels."""
    window, mask = circle_mask(src, x, y, radius)
    if window is None:
        return None, 0

    red = src.read(BAND_RED, window=window).astype("float32")
    nir = src.read(BAND_NIR, window=window).astype("float32")
    alpha = src.read(BAND_MS_ALPHA, window=window)

    denom = nir + red
    valid = mask & (alpha > 0) & (denom > 0)
    if not valid.any():
        return None, 0

    ndvi = (nir[valid] - red[valid]) / denom[valid]
    return float(ndvi.mean()), int(valid.sum())


def extract_thumbnail(src, x, y, radius, out_path, labels, font):
    """Write the annotated visible chip covering the buffer bbox.

    Draws a faint open circle of RING_DIAMETER at the survey point and a
    translucent label box in the lower left. Returns True on success.
    """
    window = from_bounds(
        x - radius, y - radius, x + radius, y + radius, src.transform
    ).round_offsets().round_lengths()
    window = window.intersection(Window(0, 0, src.width, src.height))
    if window.width <= 0 or window.height <= 0:
        return False

    rgb = src.read((1, 2, 3), window=window)
    img = Image.fromarray(np.moveaxis(rgb, 0, -1)).convert("RGBA")
    if SUPERSAMPLE > 1:
        img = img.resize(
            (img.width * SUPERSAMPLE, img.height * SUPERSAMPLE), Image.LANCZOS
        )
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # Point position within the window, in pixels (window rounding means this
    # is not exactly the image center)
    row, col = src.index(x, y)
    cx = (col - window.col_off + 0.5) * SUPERSAMPLE
    cy = (row - window.row_off + 0.5) * SUPERSAMPLE
    ring_px = (RING_DIAMETER / 2.0) / src.res[0] * SUPERSAMPLE
    draw.ellipse(
        [cx - ring_px, cy - ring_px, cx + ring_px, cy + ring_px],
        outline=RING_COLOR,
        width=RING_WIDTH * SUPERSAMPLE,
    )

    # Translucent label box in the lower left
    pad = LABEL_PAD * SUPERSAMPLE
    spacing = LABEL_LINE_SPACING * SUPERSAMPLE
    text = "\n".join(labels)
    left, top, right, bottom = draw.multiline_textbbox(
        (0, 0), text, font=font, spacing=spacing
    )
    box_w = right - left + 2 * pad
    box_h = bottom - top + 2 * pad
    box_x = pad
    box_y = img.height - pad - box_h
    draw.rectangle([box_x, box_y, box_x + box_w, box_y + box_h], fill=LABEL_BG)
    draw.multiline_text(
        (box_x + pad - left, box_y + pad - top),
        text,
        font=font,
        fill=LABEL_FG,
        spacing=spacing,
    )

    Image.alpha_composite(img, overlay).convert("RGB").save(out_path)
    return True


def extract_shadow_index(src, x, y, radius):
    """Mean Rikimaru shadow index over a circular buffer of the visible ortho."""
    window, mask = circle_mask(src, x, y, radius)
    if window is None:
        return None

    rgb = src.read((1, 2, 3), window=window).astype("float32") / 255.0
    alpha = src.read(BAND_RGB_ALPHA, window=window)

    valid = mask & (alpha > 0)
    if not valid.any():
        return None

    si = np.cbrt(
        (1.0 - rgb[0][valid]) * (1.0 - rgb[1][valid]) * (1.0 - rgb[2][valid])
    )
    return float(si.mean())


_WORKER = {}


def init_worker():
    """Open per-process GDAL handles; datasets can't be shared across forks."""
    env = rasterio.Env(**GDAL_OPTS)
    env.__enter__()
    _WORKER["env"] = env
    _WORKER["ms"] = rasterio.open(INPUT_MS)
    _WORKER["rgb"] = rasterio.open(INPUT_RGB)
    _WORKER["font"] = load_font(LABEL_FONT_SIZE * SUPERSAMPLE)


def process_point(task):
    """Extract every value for one point and write its annotated thumbnail."""
    idx, seedling_id, x, y, gps_rms, tilt = task
    ms, rgb, font = _WORKER["ms"], _WORKER["rgb"], _WORKER["font"]

    record = {"idx": idx}
    for field, radius in NDVI_RADII.items():
        value, n_px = extract_ndvi(ms, x, y, radius)
        record[field] = value
        record[f"n_px_{field.split('_')[-1]}"] = n_px
    record["shadow_idx"] = extract_shadow_index(rgb, x, y, THUMB_RADIUS)

    labels = [
        seedling_id,
        f"NDVI 0.1 m:  {fmt(record['ndvi_0o1'])}",
        f"NDVI 0.25 m: {fmt(record['ndvi_0o25'])}",
        f"Shadow Index: {fmt(record['shadow_idx'])}",
        f"GPS RMS: {fmt_quality(gps_rms, 'm')}",
        f"Tilt: {fmt_quality(tilt, 'deg', decimals=1)}",
    ]
    thumb_path = OUTPUT_THUMBS / f"{seedling_id}.png"
    written = extract_thumbnail(rgb, x, y, THUMB_RADIUS, thumb_path, labels, font)
    record["thumbnail"] = str(thumb_path.relative_to(BASE_DIR)) if written else None

    return record


def main():
    """Main processing function."""
    print("=" * 60)
    print("NDVI and Shadow Index Extraction")
    print("Miller Creek Collaring - 2025-10-21 survey")
    print("=" * 60)

    gdf = load_points()
    OUTPUT_THUMBS.mkdir(parents=True, exist_ok=True)
    OUTPUT_GEOJSON.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.Env(**GDAL_OPTS):
        print("\nOpening orthomosaics over /vsicurl...")
        with rasterio.open(INPUT_MS) as ms, rasterio.open(INPUT_RGB) as rgb:
            print(f"  MS:  {ms.count} bands, {ms.res[0]:.4f} m GSD, {ms.crs}")
            print(f"  RGB: {rgb.count} bands, {rgb.res[0]:.4f} m GSD, {rgb.crs}")

            if ms.crs != rgb.crs:
                raise ValueError("MS and RGB orthomosaics have different CRS")

            # Project points into the ortho CRS for metric buffering
            gdf_ortho = gdf.to_crs(ms.crs)

    tasks = [
        (
            int(row["idx"]),
            str(row["Name"]),
            row.geometry.x,
            row.geometry.y,
            row["gps_rms_m"],
            row["tilt_deg"],
        )
        for _, row in gdf_ortho.iterrows()
    ]

    print(f"\nExtracting {len(tasks)} points on {N_WORKERS} workers...")
    records = []
    with ProcessPoolExecutor(
        max_workers=N_WORKERS, initializer=init_worker
    ) as pool:
        for i, record in enumerate(pool.map(process_point, tasks), start=1):
            records.append(record)
            if i % 20 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)}")

    out = gdf[[
        "idx", "Name", "Code description", "Description",
        "gps_rms_m", "tilt_deg", "geometry",
    ]].merge(pd.DataFrame(records).sort_values("idx"), on="idx")
    out = out.rename(columns={"Code description": "code_desc", "Description": "note"})

    # Keep the requested fields up front
    cols = [
        "idx", "Name", "code_desc", "note",
        "ndvi_0o25", "ndvi_0o1", "shadow_idx",
        "gps_rms_m", "tilt_deg",
        "n_px_0o25", "n_px_0o1", "thumbnail", "geometry",
    ]
    out = gpd.GeoDataFrame(out[cols], geometry="geometry", crs=gdf.crs)

    out.to_file(OUTPUT_GEOJSON, driver="GeoJSON")

    print("\nSummary:")
    for field in ["ndvi_0o25", "ndvi_0o1", "shadow_idx"]:
        values = out[field].dropna()
        print(
            f"  {field:11s} n={len(values):3d}  "
            f"min={values.min():+.3f}  mean={values.mean():+.3f}  max={values.max():+.3f}"
        )
    for field, unit in [("gps_rms_m", "m"), ("tilt_deg", "deg")]:
        values = out[field].dropna()
        absent = len(out) - len(values)
        flag = f"  [{absent} x {MISSING_LABEL}]" if absent else ""
        print(
            f"  {field:11s} n={len(values):3d}  "
            f"min={values.min():.3f}  mean={values.mean():.3f}  "
            f"max={values.max():.3f} {unit}{flag}"
        )

    missing = int(out["ndvi_0o25"].isna().sum())
    if missing:
        print(f"  {missing} point(s) fell outside the orthomosaic footprint")

    print("\n" + "=" * 60)
    print("Processing complete!")
    print(f"Thumbnails: {OUTPUT_THUMBS}")
    print(f"Output saved to: {OUTPUT_GEOJSON}")
    print("=" * 60)


if __name__ == "__main__":
    main()
