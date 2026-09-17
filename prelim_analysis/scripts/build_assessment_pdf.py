#!/usr/bin/env python3
"""
Join the annotated seedling thumbnails into a single assessment PDF.

Builds a cover page (sources, method, caveats, summary statistics) followed by
contact-sheet pages of every annotated chip, ordered by seedling ID.

Pages are rasterized at DPI; every layout constant and font size is expressed at
150 dpi and scaled through S(), so raising DPI sharpens the text without moving
anything on the page.

Inputs:
    - prelim_analysis/results/ndvi_shadows.geojson: extracted values
    - data/raster/seedling_thumbnails/2025/{seedling_id}.png: annotated chips

Outputs:
    - prelim_analysis/initial_uncalibrated_ndvi_assessment.pdf
"""

import geopandas as gpd
import sys
from pathlib import Path
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).parent))
from ndvi_extract import load_font

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_GEOJSON = BASE_DIR / "prelim_analysis" / "results" / "ndvi_shadows.geojson"
OUTPUT_PDF = (
    BASE_DIR / "prelim_analysis" / "initial_uncalibrated_ndvi_assessment.pdf"
)

# US Letter. Layout below is authored at 150 dpi and scaled by S().
DPI = 300
BASE_DPI = 150
SCALE = DPI / BASE_DPI
PAGE_W, PAGE_H = int(8.5 * DPI), int(11 * DPI)

MARGIN = 40
GAP = 22
HEADER_H = 30
COLS, ROWS = 2, 3

SURVEY_DATE = "2025-10-21"
MS_GSD_CM = 2.66
RGB_GSD_CM = 1.28

INK = (25, 25, 25)
MUTED = (110, 110, 110)
RULE = (200, 200, 200)

MISSING_LABEL = "!!!NONE!!!"


def S(value):
    """Scale a 150-dpi layout value to the output resolution."""
    return int(round(value * SCALE))


def new_page():
    """Blank letter page."""
    return Image.new("RGB", (PAGE_W, PAGE_H), "white")


def build_cover(gdf, fonts):
    """Cover page: sources, method, caveats, and summary statistics."""
    page = new_page()
    draw = ImageDraw.Draw(page)
    y = S(MARGIN + 20)

    draw.text((S(MARGIN), y), "Miller Creek Collaring", font=fonts["title"], fill=INK)
    y += S(46)
    draw.text(
        (S(MARGIN), y),
        "Initial uncalibrated NDVI and shadow assessment",
        font=fonts["subtitle"],
        fill=INK,
    )
    y += S(34)
    draw.text(
        (S(MARGIN), y),
        f"Survey flown {SURVEY_DATE}  |  {len(gdf)} surveyed points",
        font=fonts["body"],
        fill=MUTED,
    )
    y += S(40)
    draw.line([(S(MARGIN), y), (PAGE_W - S(MARGIN), y)], fill=RULE, width=S(2))
    y += S(26)

    sections = [
        ("Source imagery", [
            "DroneDeploy orthomosaics, 50 m AGL, GCP-controlled, EPSG:6514",
            f"Multispectral: 5 bands + alpha, {MS_GSD_CM} cm GSD (R, G, NIR, RE, unused, alpha)",
            f"Visible: RGB + alpha, {RGB_GSD_CM} cm GSD",
            "Survey points: Emlid RS3 RTK, all fixed solutions, EPSG:4326",
            "GPS RMS below is the Lateral (horizontal) RMS from the Emlid export;",
            "tilt is the pole tilt angle. Absent values are marked !!!NONE!!!",
        ]),
        ("Method", [
            "NDVI = (NIR - Red) / (NIR + Red), band 3 and band 1 of the MS ortho",
            "Averaged per pixel over circular buffers of 0.10 m and 0.25 m radius",
            "(44 and ~277 contributing pixels respectively; alpha = 0 excluded)",
            "Shadow index after Rikimaru et al. (2002): SI = [(1-R)(1-G)(1-B)]^(1/3)",
            "on 0-1 scaled visible bands, averaged over a 3 m radius buffer.",
            "Higher SI means darker / more shadowed.",
        ]),
        ("Caveat", [
            "The MS ortho carries isCalibrated = False and 8-bit bands, so NDVI here",
            "is a relative index on raw DN, not reflectance NDVI. Values are",
            "comparable within this survey only, not across dates or sensors.",
        ]),
    ]
    for heading, lines in sections:
        draw.text((S(MARGIN), y), heading, font=fonts["heading"], fill=INK)
        y += S(26)
        for line in lines:
            draw.text((S(MARGIN + 16), y), line, font=fonts["body"], fill=INK)
            y += S(21)
        y += S(16)

    # Summary statistics
    draw.text((S(MARGIN), y), "Summary statistics", font=fonts["heading"], fill=INK)
    y += S(28)
    cols = [S(MARGIN + off) for off in (16, 200, 300, 400, 500, 600)]
    headers = ["field", "n", "min", "mean", "median", "max"]
    for x, label in zip(cols, headers):
        draw.text((x, y), label, font=fonts["mono_bold"], fill=MUTED)
    y += S(22)
    draw.line([(S(MARGIN + 16), y), (PAGE_W - S(MARGIN), y)], fill=RULE, width=S(1))
    y += S(8)

    stat_fields = [
        ("ndvi_0o1", "{:+.3f}"),
        ("ndvi_0o25", "{:+.3f}"),
        ("shadow_idx", "{:+.3f}"),
        ("gps_rms_m", "{:.3f}"),
        ("tilt_deg", "{:.1f}"),
    ]
    n_absent = {}
    for field, spec in stat_fields:
        values = gdf[field].dropna()
        n_absent[field] = len(gdf) - len(values)
        cells = [field, f"{len(values)}"] + [
            spec.format(v)
            for v in (values.min(), values.mean(), values.median(), values.max())
        ]
        for x, cell in zip(cols, cells):
            draw.text((x, y), cell, font=fonts["mono"], fill=INK)
        y += S(24)

    y += S(8)
    for field, absent in n_absent.items():
        if absent:
            draw.text(
                (S(MARGIN + 16), y),
                f"{absent} of {len(gdf)} points have no {field} recorded; "
                f"marked {MISSING_LABEL} on the chips.",
                font=fonts["body"],
                fill=INK,
            )
            y += S(21)

    draw.text(
        (S(MARGIN), PAGE_H - S(MARGIN + 18)),
        "Generated by prelim_analysis/scripts/build_assessment_pdf.py",
        font=fonts["small"],
        fill=MUTED,
    )
    return page


def build_sheets(gdf, fonts):
    """Contact-sheet pages, COLS x ROWS annotated chips per page."""
    records = [r for _, r in gdf.iterrows() if r["thumbnail"]]
    per_page = COLS * ROWS
    n_pages = (len(records) + per_page - 1) // per_page

    margin, gap, header_h = S(MARGIN), S(GAP), S(HEADER_H)
    cell_w = (PAGE_W - 2 * margin - (COLS - 1) * gap) // COLS
    cell_h = (PAGE_H - 2 * margin - header_h - (ROWS - 1) * gap) // ROWS
    size = min(cell_w, cell_h)
    grid_w = COLS * size + (COLS - 1) * gap
    x0 = (PAGE_W - grid_w) // 2

    pages = []
    for p in range(n_pages):
        page = new_page()
        draw = ImageDraw.Draw(page)
        draw.text(
            (margin, margin - S(12)),
            f"Miller Creek collaring  |  uncalibrated NDVI assessment  |  "
            f"{SURVEY_DATE}  |  sheet {p + 1} of {n_pages}",
            font=fonts["small"],
            fill=MUTED,
        )

        chunk = records[p * per_page:(p + 1) * per_page]
        for k, record in enumerate(chunk):
            chip = Image.open(BASE_DIR / record["thumbnail"]).convert("RGB")
            chip = chip.resize((size, size), Image.LANCZOS)
            x = x0 + (k % COLS) * (size + gap)
            y = margin + header_h + (k // COLS) * (size + gap)
            page.paste(chip, (x, y))
            draw.rectangle(
                [x, y, x + size - 1, y + size - 1], outline=RULE, width=S(1)
            )

        pages.append(page)

    return size, pages


def main():
    """Main processing function."""
    print("=" * 60)
    print("Building NDVI Assessment PDF")
    print("Miller Creek Collaring")
    print("=" * 60)

    print("\nLoading extracted values...")
    gdf = gpd.read_file(INPUT_GEOJSON).sort_values("Name").reset_index(drop=True)
    print(f"  {len(gdf)} points loaded")

    fonts = {
        "title": load_font(S(34)),
        "subtitle": load_font(S(22)),
        "heading": load_font(S(18)),
        "body": load_font(S(14)),
        "small": load_font(S(12)),
        "mono": load_font(S(14)),
        "mono_bold": load_font(S(14)),
    }

    print(f"\nRendering pages at {DPI} dpi ({PAGE_W}x{PAGE_H} px)...")
    cover = build_cover(gdf, fonts)
    cell_px, sheets = build_sheets(gdf, fonts)
    chip_px = Image.open(BASE_DIR / gdf.iloc[0]["thumbnail"]).width
    print(f"  1 cover + {len(sheets)} contact sheets")
    print(f"  chips {chip_px} px rendered into {cell_px} px cells")

    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    cover.save(
        OUTPUT_PDF,
        save_all=True,
        append_images=sheets,
        resolution=DPI,
        title="Miller Creek initial uncalibrated NDVI assessment",
    )

    size_mb = OUTPUT_PDF.stat().st_size / 1e6

    print("\n" + "=" * 60)
    print("Processing complete!")
    print(f"Output saved to: {OUTPUT_PDF} ({size_mb:.1f} MB, {1 + len(sheets)} pages)")
    print("=" * 60)


if __name__ == "__main__":
    main()
