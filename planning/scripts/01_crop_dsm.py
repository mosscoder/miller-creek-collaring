#!/usr/bin/env python3
"""
Crop DSM raster to buffered mission boundary and convert to ellipsoidal height.

Inputs:
    - data/vector/miller_collaring_mission.kml: Mission boundary (30m buffer)
    - data/raster/lidar_products/46113g8/MISSOULA_2021_Phase3_1/DSM/46113g8_DSM.tif: Orthometric DSM
    - data/raster/g2018u0.bin: Geoid model

Outputs:
    - data/raster/miller_terrain_follow.tif: Cropped ellipsoidal DSM in EPSG:4326
"""

import geopandas as gpd
import rasterio
import subprocess
import json
from pathlib import Path

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_KML = BASE_DIR / "data" / "vector" / "miller_collaring_mission.kml"
INPUT_DSM = BASE_DIR / "data" / "raster" / "lidar_products" / "46113g8" / "MISSOULA_2021_Phase3_1" / "DSM" / "46113g8_DSM.tif"
INPUT_GEOID = BASE_DIR / "data" / "raster" / "g2018u0.bin"
OUTPUT_DSM = BASE_DIR / "data" / "raster" / "miller_terrain_follow.tif"

# Intermediate files
TEMP_DIR = BASE_DIR / "data" / "raster"
CROPPED_DSM = TEMP_DIR / "temp_cropped_dsm.tif"
WARPED_GEOID = TEMP_DIR / "temp_warped_geoid.tif"
ELLIPSOIDAL_DSM_NATIVE = TEMP_DIR / "temp_ellipsoidal_native.tif"

def main():
    """Main processing function."""
    print("=" * 60)
    print("Cropping and Converting DSM to Ellipsoidal Height")
    print("Miller Creek Collaring")
    print("=" * 60)

    print("\nLoading mission boundary...")
    # Read mission boundary KML
    gdf = gpd.read_file(INPUT_KML)

    # Ensure EPSG:4326
    if gdf.crs is None:
        gdf.set_crs("EPSG:4326", inplace=True)
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")

    print("  Mission boundary loaded")

    print("\nApplying 50m buffer...")
    # Project to UTM Zone 12N for accurate metric buffering
    gdf_utm = gdf.to_crs("EPSG:32612")

    # Buffer by 50 meters
    buffered = gdf_utm.buffer(50)
    gdf_buffered = gpd.GeoDataFrame(geometry=buffered, crs="EPSG:32612")

    print(f"  Buffered by 50m in UTM Zone 12N")

    # Get DSM properties
    print("\nExtracting DSM parameters...")
    with rasterio.open(INPUT_DSM) as src:
        dsm_crs = src.crs
        print(f"  DSM CRS: {dsm_crs}")

        # Reproject buffer to DSM CRS to get bounds
        gdf_dsm_crs = gdf_buffered.to_crs(dsm_crs)
        bounds = gdf_dsm_crs.total_bounds
        xmin, ymin, xmax, ymax = bounds

        # Get resolution
        xres, yres = src.res

    print(f"  Bounds: {xmin} {ymin} {xmax} {ymax}")
    print(f"  Resolution: {xres} {yres}")

    # Step 1: Crop DSM to buffered boundary
    print("\nCropping DSM to buffered boundary...")
    subprocess.run([
        'gdalwarp',
        '-r', 'bilinear',
        '-te', str(xmin), str(ymin), str(xmax), str(ymax),
        '-tr', str(xres), str(abs(yres)),
        '-crop_to_cutline',
        str(INPUT_DSM), str(CROPPED_DSM)
    ], check=True)
    print(f"  Cropped DSM saved: {CROPPED_DSM}")

    # Step 2: Warp geoid to match cropped DSM
    print("\nWarping geoid to match cropped DSM...")
    subprocess.run([
        'gdalwarp',
        '-r', 'bilinear',
        '-t_srs', str(dsm_crs),
        '-te', str(xmin), str(ymin), str(xmax), str(ymax),
        '-tr', str(xres), str(abs(yres)),
        '-crop_to_cutline',
        str(INPUT_GEOID), str(WARPED_GEOID)
    ], check=True)
    print(f"  Warped geoid saved: {WARPED_GEOID}")

    # Step 3: Convert to ellipsoidal height (Ellipsoidal = Orthometric + Geoid)
    print("\nConverting to ellipsoidal height...")
    subprocess.run([
        'gdal_calc.py',
        '--calc=A+B',
        '--format=GTiff',
        '--type=Float32',
        '-A', str(CROPPED_DSM), '--A_band=1',
        '-B', str(WARPED_GEOID), '--B_band=1',
        '--outfile=' + str(ELLIPSOIDAL_DSM_NATIVE),
        '--co=COMPRESS=DEFLATE',
        '--co=PREDICTOR=2',
        '--co=TILED=YES',
        '--overwrite'
    ], check=True)
    print(f"  Ellipsoidal DSM (native CRS) saved: {ELLIPSOIDAL_DSM_NATIVE}")

    # Step 4: Reproject to EPSG:4326
    print("\nReprojecting to EPSG:4326...")
    subprocess.run([
        'gdalwarp',
        '-r', 'bilinear',
        '-t_srs', 'EPSG:4326',
        '-co', 'COMPRESS=DEFLATE',
        '-co', 'PREDICTOR=2',
        '-co', 'TILED=YES',
        str(ELLIPSOIDAL_DSM_NATIVE), str(OUTPUT_DSM)
    ], check=True)
    print(f"  Final output saved: {OUTPUT_DSM}")

    # Step 5: Clean up intermediate files
    print("\nCleaning up intermediate files...")
    CROPPED_DSM.unlink()
    WARPED_GEOID.unlink()
    ELLIPSOIDAL_DSM_NATIVE.unlink()
    print("  Cleanup complete")

    print("\n" + "=" * 60)
    print("Processing complete!")
    print(f"Output saved to: {OUTPUT_DSM}")
    print("=" * 60)

if __name__ == "__main__":
    main()
