## Project Overview
This project contains workflows and spatial data for surveying regions in Miller Creek, Montana, focused on capturing tree images before collaring and subsequent death.

## Action items
- Get Beau's areas of interest, do they align with Mary-Ellyn's AOI
- Day 1 - Bagging GCPs with LoRA
- Day 2 - Flying after GCPs acquired

## Key files:
- `data/vector/miller_collaring_mission.kml` -  Buffered mission boundary for DJI Pilot 2
- `data/vector/gcps_init.csv` - GCPs formatted for Emlid Flow (EPSG:4326)
- `data/raster/miller_terrain_follow.tif` - ellipsoidal DSM in EPSG:4326 for DJI Pilot 2 terrain following

## Mission Map:
![Mission Map](planning/miller_mission_map.png)

## Flight Specs:
- AGL: 50m
- Sensor: RGB + MS
- Shutter speed: 1/1000s
- Lap: 80/70
- Flight Speed: Auto
- Terrain following: Enabled
- GCP count: 6

## Drone Deploy Exports:
- RGB orthomosaic
- 5-band MS orthomosaic
- .las Point Cloud
