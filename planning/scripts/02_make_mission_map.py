#!/usr/bin/env python3
"""
Create static mission map with intensity, CHM, and vector overlays.

Inputs:
    - data/raster/lidar_products/46113g8/MISSOULA_2021_Phase3_1/Intensity/46113g8_Intensity.tif: LiDAR intensity
    - data/raster/lidar_products/46113g8/MISSOULA_2021_Phase3_1/CHM/46113g8_CHM.tif: Canopy Height Model
    - data/vector/init_soil_team_aoi.kml: Initial AOI
    - data/vector/miller_collaring_mission.kml: Buffered AOI
    - data/vector/gcps_init.geojson: GCP points

Outputs:
    - planning/miller_mission_map.png: Static map
"""

import geopandas as gpd
import rasterio
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from rasterio.plot import show
from pathlib import Path
import numpy as np

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_HILLSHADE = BASE_DIR / "data" / "raster" / "lidar_products" / "46113g8" / "MISSOULA_2021_Phase3_1" / "Intensity" / "46113g8_Intensity.tif"
INPUT_CHM = BASE_DIR / "data" / "raster" / "lidar_products" / "46113g8" / "MISSOULA_2021_Phase3_1" / "CHM" / "46113g8_CHM.tif"
INPUT_AOI_INIT = BASE_DIR / "data" / "vector" / "init_soil_team_aoi.kml"
INPUT_AOI_BUFFERED = BASE_DIR / "data" / "vector" / "miller_collaring_mission.kml"
INPUT_GCPS = BASE_DIR / "data" / "vector" / "gcps_init.geojson"

CROPPED_HILLSHADE = BASE_DIR / "data" / "raster" / "miller_terrain_follow_hillshade.tif"
HILLSHADE_6514 = BASE_DIR / "data" / "raster" / "temp_hillshade_6514.tif"
CROPPED_CHM = BASE_DIR / "data" / "raster" / "temp_cropped_chm.tif"
CHM_6514 = BASE_DIR / "data" / "raster" / "temp_chm_6514.tif"
MAP_OUTPUT = BASE_DIR / "planning" / "miller_mission_map.png"

def main():
    """Main processing function."""
    print("=" * 60)
    print("Creating Mission Map")
    print("Miller Creek Collaring")
    print("=" * 60)

    # Step 1: Load and reproject vector data
    print("\nLoading vector data...")

    # Load initial AOI
    aoi_init = gpd.read_file(INPUT_AOI_INIT)
    if aoi_init.crs is None:
        aoi_init.set_crs("EPSG:4326", inplace=True)
    aoi_init = aoi_init.to_crs("EPSG:6514")
    print(f"  Initial AOI loaded and reprojected")

    # Load buffered AOI
    aoi_buffered = gpd.read_file(INPUT_AOI_BUFFERED)
    if aoi_buffered.crs is None:
        aoi_buffered.set_crs("EPSG:4326", inplace=True)
    aoi_buffered = aoi_buffered.to_crs("EPSG:6514")
    print(f"  Buffered AOI loaded and reprojected")

    # Load GCPs
    gcps = gpd.read_file(INPUT_GCPS)
    if gcps.crs is None:
        gcps.set_crs("EPSG:4326", inplace=True)
    gcps = gcps.to_crs("EPSG:6514")
    print(f"  GCPs loaded and reprojected ({len(gcps)} points)")

    # Step 2: Crop and reproject intensity
    print("\nProcessing intensity...")

    # Get intensity CRS and get buffered boundary in that CRS
    with rasterio.open(INPUT_HILLSHADE) as src:
        intensity_crs = src.crs
        xres, yres = src.res

    # Reproject buffered boundary to intensity CRS to get crop bounds
    aoi_buffered_intensity_crs = aoi_buffered.to_crs(intensity_crs)
    bounds = aoi_buffered_intensity_crs.total_bounds
    xmin, ymin, xmax, ymax = bounds

    # Extend bounds by 20 meters for plot footprint
    xmin -= 20
    ymin -= 20
    xmax += 20
    ymax += 20

    print(f"  Intensity CRS: {intensity_crs}")
    print(f"  Cropping to buffered boundary (+20m extension)...")

    # Crop intensity to buffered boundary
    subprocess.run([
        'gdalwarp',
        '-r', 'bilinear',
        '-te', str(xmin), str(ymin), str(xmax), str(ymax),
        '-tr', str(xres), str(abs(yres)),
        str(INPUT_HILLSHADE), str(CROPPED_HILLSHADE)
    ], check=True)

    # Reproject to EPSG:6514
    subprocess.run([
        'gdalwarp',
        '-t_srs', 'EPSG:6514',
        '-r', 'bilinear',
        str(CROPPED_HILLSHADE),
        str(HILLSHADE_6514)
    ], check=True)

    # Step 3: Crop and reproject CHM
    print("\nProcessing CHM...")
    print(f"  Cropping to buffered boundary (+20m extension)...")

    # Crop CHM to buffered boundary (same bounds as intensity)
    subprocess.run([
        'gdalwarp',
        '-r', 'bilinear',
        '-te', str(xmin), str(ymin), str(xmax), str(ymax),
        '-tr', str(xres), str(abs(yres)),
        str(INPUT_CHM), str(CROPPED_CHM)
    ], check=True)

    # Reproject to EPSG:6514
    subprocess.run([
        'gdalwarp',
        '-t_srs', 'EPSG:6514',
        '-r', 'bilinear',
        str(CROPPED_CHM),
        str(CHM_6514)
    ], check=True)

    # Step 4: Create map
    print("\nCreating map...")
    fig, ax = plt.subplots(figsize=(12, 14))

    # Plot intensity
    with rasterio.open(HILLSHADE_6514) as src:
        show(src, ax=ax, cmap='gray', alpha=1.0)
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

    # Plot CHM (viridis colormap, 0.8 transparency, mask < 1.0m)
    with rasterio.open(CHM_6514) as src:
        chm_data = src.read(1)
        chm_nodata = src.nodata
        chm_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

        # Mask CHM values < 1.0m and nodata
        chm_masked = np.ma.masked_where(chm_data < 1.0, chm_data)
        if chm_nodata is not None:
            chm_masked = np.ma.masked_where(chm_data == chm_nodata, chm_masked)

        # Plot masked CHM
        ax.imshow(chm_masked, cmap='viridis', alpha=0.5,
                 extent=chm_extent, origin='upper', interpolation='bilinear')

    # Plot initial AOI (red)
    aoi_init.boundary.plot(ax=ax, color='red', linewidth=2, label='Initial AOI')

    # Plot buffered AOI (cyan)
    aoi_buffered.boundary.plot(ax=ax, color='orange', linewidth=2, label='Buffered Mission Boundary')

    # Plot GCPs (red cross on white dot)
    gcps_coords = [(point.x, point.y) for point in gcps.geometry]
    x_coords = [coord[0] for coord in gcps_coords]
    y_coords = [coord[1] for coord in gcps_coords]

    # White dots
    ax.scatter(x_coords, y_coords, c='white', s=100, zorder=5, edgecolors='black', linewidth=1)
    # Red crosses
    ax.scatter(x_coords, y_coords, c='red', s=50, marker='x', zorder=6, linewidth=2)

    # Add grid
    ax.grid(True, linestyle='--', alpha=0.5, color='gray')
    ax.set_xlabel('Easting (m)', fontsize=12)
    ax.set_ylabel('Northing (m)', fontsize=12)
    ax.set_title('Miller Creek Collaring Mission Map', fontsize=14, fontweight='bold')

    # Create legend
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, label='Initial AOI'),
        Line2D([0], [0], color='orange', linewidth=2, label='Buffered Mission Boundary'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='white',
               markeredgecolor='black', markersize=10, label='GCPs',
               markeredgewidth=1, linestyle='None')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    # Set aspect ratio to equal
    ax.set_aspect('equal')

    # Tight layout
    plt.tight_layout()

    # Save map
    print(f"\nSaving map to {MAP_OUTPUT}...")
    plt.savefig(MAP_OUTPUT, dpi=300, bbox_inches='tight')
    print(f"  Map saved successfully")

    # Clean up temporary files
    print("\nCleaning up temporary files...")
    CROPPED_HILLSHADE.unlink()
    HILLSHADE_6514.unlink()
    CROPPED_CHM.unlink()
    CHM_6514.unlink()
    print("  Cleanup complete")

    print("\n" + "=" * 60)
    print("Processing complete!")
    print(f"Map: {MAP_OUTPUT}")
    print("=" * 60)

if __name__ == "__main__":
    main()
