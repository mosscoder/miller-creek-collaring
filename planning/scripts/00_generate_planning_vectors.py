#!/usr/bin/env python3
"""
Generate planning vectors for Miller Creek collaring project.

Inputs:
    - data/vector/gcps_init.geojson: Initial GCP points
    - data/vector/init_soil_team_aoi.kml: Initial AOI polygon

Outputs:
    - data/vector/gcps_init.csv: GCPs formatted for Emlid Flow (EPSG:4326)
    - data/vector/miller_collaring_mission.kml: Buffered mission boundary for DJI Pilot 2
"""

import geopandas as gpd
import pandas as pd
from pathlib import Path

# Define paths
BASE_DIR = Path(__file__).parent.parent.parent
INPUT_GCP = BASE_DIR / "data" / "vector" / "gcps_init.geojson"
INPUT_AOI = BASE_DIR / "data" / "vector" / "init_soil_team_aoi.kml"
OUTPUT_GCP = BASE_DIR / "data" / "vector" / "gcps_init.csv"
OUTPUT_KML = BASE_DIR / "data" / "vector" / "miller_collaring_mission.kml"

def process_gcps():
    """Process GCP points for Emlid Flow format."""
    print("Processing GCPs...")

    # Read GeoJSON
    gdf = gpd.read_file(INPUT_GCP)

    # Ensure EPSG:4326
    if gdf.crs is None:
        gdf.set_crs("EPSG:4326", inplace=True)
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")

    # Extract coordinates
    gdf['longitude'] = gdf.geometry.x
    gdf['latitude'] = gdf.geometry.y

    # Sort: upper-left to lower-right (latitude descending, then longitude ascending)
    gdf = gdf.sort_values(['latitude', 'longitude'], ascending=[False, True])
    gdf = gdf.reset_index(drop=True)

    # Generate GCP IDs
    gdf['name'] = [f"GCP_{i+1:03d}" for i in range(len(gdf))]

    # Set altitude to 0 (will be measured in field)
    gdf['ellipsoid_height'] = 0

    # Create output dataframe
    output_df = gdf[['name', 'longitude', 'latitude', 'ellipsoid_height']]

    # Save to CSV
    output_df.to_csv(OUTPUT_GCP, index=False)
    print(f"  Saved {len(output_df)} GCPs to {OUTPUT_GCP}")
    print(f"  GCP ordering: {', '.join(output_df['name'].tolist())}")

def process_aoi():
    """Process AOI polygon with 30m buffer for DJI Pilot 2."""
    print("\nProcessing AOI...")

    # Read KML
    gdf = gpd.read_file(INPUT_AOI)

    # Ensure EPSG:4326
    if gdf.crs is None:
        gdf.set_crs("EPSG:4326", inplace=True)
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")

    # Project to UTM Zone 12N (EPSG:32612) for accurate metric buffering
    # Montana at longitude -113.94° falls in UTM Zone 12N
    gdf_utm = gdf.to_crs("EPSG:32612")

    # Buffer by 50 meters with miter join style (sharp corners)
    buffered = gdf_utm.buffer(30, join_style=2)  # join_style=2 is miter

    # Create new GeoDataFrame with buffered geometry
    gdf_buffered = gpd.GeoDataFrame(geometry=buffered, crs="EPSG:32612")

    # Transform back to WGS84
    gdf_buffered = gdf_buffered.to_crs("EPSG:4326")

    # Generate simple KML for DJI Pilot 2
    write_simple_kml(gdf_buffered, OUTPUT_KML)
    print(f"  Saved buffered mission boundary to {OUTPUT_KML}")

def write_simple_kml(gdf, output_path):
    """
    Write a minimal KML file compatible with DJI Pilot 2.

    DJI Pilot 2 requires simple KML structure without complex styling or nesting.
    """
    # Get the first geometry (should only be one polygon)
    geom = gdf.geometry.iloc[0]

    # Extract exterior coordinates
    coords = geom.exterior.coords

    # Format coordinates as lon,lat,0 (KML format)
    coord_str = " ".join([f"{lon},{lat},0" for lon, lat in coords])

    # Create minimal KML structure
    kml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
<name>Miller Creek Collaring Mission</name>
<Placemark>
<name>Mission Boundary</name>
<Polygon>
<outerBoundaryIs>
<LinearRing>
<coordinates>
{coord_str}
</coordinates>
</LinearRing>
</outerBoundaryIs>
</Polygon>
</Placemark>
</Document>
</kml>
"""

    # Write to file
    with open(output_path, 'w') as f:
        f.write(kml_content)

def main():
    """Main processing function."""
    print("=" * 60)
    print("Generating Planning Vectors for Miller Creek Collaring")
    print("=" * 60)

    # Process GCPs
    process_gcps()

    # Process AOI
    process_aoi()

    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
