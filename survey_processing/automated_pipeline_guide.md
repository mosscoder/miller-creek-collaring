# Processing the 2026-09-08 collaring flight with the automated pipeline

The `aerial` module takes a flight's raw frames from the survey bucket, builds the DroneDeploy uploads
in the cloud (the multispectral frames get DJI's radiometric calibration on the way), creates the maps in
DroneDeploy, waits for them to process, lands the exports back in the bucket, and then registers the
multispectral map to the visible map on a GPU. Nothing passes through your laptop except this config
file and a few commands. The module lives at https://github.com/mosscoder/mpg-aerial-pipeline
(`docs/operations.md` there explains every moving part).

What you will end up with, under `gs://mpg-aerial-survey/surveys/miller_collaring/260908/processing/drone_deploy/`:

```
miller_collaring-260908-visible.tif             RGBA orthomosaic, GCP-rectified, EPSG 6514
miller_collaring-260908-pointcloud.las
miller_collaring-260908-multispectral.tif       Red, Green, NIR, RedEdge + alpha, registered to the visible map
miller_collaring-260908-*.manifest.json         every frame uploaded, md5s, the calibration terms
miller_collaring-260908-*.run.json              the state of each map, step by step
registration_results/miller_collaring-260908/   quality.tif, qa.json, quicklook.png, report.html
survey_config.toml                              the working copy of the file you write below
```

## 1. Set up once

Python 3.10 or newer, in whatever environment you like:

```bash
pip install git+https://github.com/mosscoder/mpg-aerial-pipeline@v0.4.0
```

Then, in this folder (`survey_processing/`), a `.env` file with two lines. Both keys stay on your machine;
the file is gitignored.

```
GOOGLE_APPLICATION_CREDENTIALS=/path/to/your-service-account-key.json
DRONE_DEPLOY_API_KEY=<the DroneDeploy API key>
```

The module refuses to run as a personal gcloud login, so the first line must point at a service-account
key. You do not need `gcloud` installed for any of this.

## 2. Put the ground control file in the bucket

DroneDeploy takes ground control as a CSV with the header `Label,Latitude,Longitude,Elevation`, WGS84
degrees and ellipsoidal metres. Save the six points that way and upload the file to the path the config
names below (any name is fine as long as the two match):

```
gs://mpg-aerial-survey/surveys/miller_collaring/260908/data_collection/ground_truth/miller_collaring-260908-gcps-dronedeploy.csv
```

The preflight in step 4 checks that the file is reachable before anything is created.

## 3. Write the config

Save this as `survey_config.toml` in this folder. It is already checked against the bucket: the four
mission folders hold 6,620 multispectral and 1,655 visible frames, and the DroneDeploy path resolves to
the existing Miller Creek project. Only the GCP file line needs to match what you uploaded.

```toml
[survey]
id            = "miller_collaring_260908"        # unique name: labels the run records and the log lines
root          = "gs://mpg-aerial-survey/surveys/miller_collaring/260908"   # where the raw flight lives; read only
raw           = "data_collection/DCIM"           # the folder under root that holds the DJI flight folders
stem          = "miller_collaring-{yymmdd}"      # name of everything this flight produces; {yymmdd} comes from the map's date
outputs       = "gs://mpg-aerial-survey/surveys/miller_collaring/260908/processing/drone_deploy"   # where products and records land
multispectral = true                             # build the four-band multispectral map (DJI-calibrated frames)
visible       = true                             # build the RGB orthomosaic and the point cloud

[dronedeploy]
path = "Side Projects / Miller Creek"            # folder / project in DroneDeploy; the plans are named <stem>-visible and <stem>-multispectral

[visible]
gcps = "data_collection/ground_truth/miller_collaring-260908-gcps-dronedeploy.csv"   # the GCP CSV, relative to root; sent with the visible upload

[maps.m260908]                                   # one block per map; m260908 is the id you use on the command line
date           = "2026-09-08"                    # fills {yymmdd} above
flight_folders = ["DJI_202609081316_021_millercollaringmission",   # the four battery folders of the mission, in order
                  "DJI_202609081342_022_millercollaringmission",
                  "DJI_202609081408_024_millercollaringmission",
                  "DJI_202609081437_025_millercollaringmission"]
notes = "folders 020, 023 and 026 are single test captures and are left out; panel_captures holds the panel shots"
```

Not listed on purpose: the three tiny folders and `panel_captures`. The module only reads what the map
names. Everything else (frame patterns, export layers, projection, the calibration rule) is the module's
default; `docs/survey_config.md` in the module repo lists every key if you ever need to change one.

## 4. Run it

```bash
aerial submit --dry-run     # preflight: counts the frames, resolves the DroneDeploy path, checks the GCP file; creates nothing
aerial submit               # builds both uploads in the cloud (about 10 min), creates the two plans, posts the transfers
aerial status               # one line per map and data type, from the run records in the bucket
```

`submit` also adds this config to the watch list of the poll job, which runs every 15 minutes in Cloud
Run and drives everything from there: it waits for DroneDeploy to fetch the uploads, follows the two maps
through processing, requests the exports, lands them in the bucket, and finally runs the registration on
the GPU job. You can close the laptop after `submit`.

One step needs a person. Ground control makes DroneDeploy tag the targets in the images and then wait
for someone to review the tags: open the visible plan in the DroneDeploy app, check the tags, and press
Continue to Processing. `aerial status` says `gcp=pending` while it waits and flags the map after two
hours. The multispectral map has no ground control and processes on its own.

## 5. Read the result

Run `aerial status` whenever you like. A map goes `submitted` → `transferred` → `processing` →
`processed` → `exporting` → `exported` → `done`; the multispectral line shows `steps: register` at the
end, and its `registration_results/` folder holds the evidence: `report.html` is self-contained, so
download it and open it in a browser. `quality.tif` is one small raster on a 64-pixel cell grid with the
measured shifts, match coverage and the agreement with the visible map before and after; `qa.json` has
the same numbers plus timing.

Expect about a day end to end, most of it DroneDeploy's queue, and roughly $10 of cost, almost all of it
DroneDeploy pulling the 80 GB of uploads out of the bucket.

## If something looks wrong

- `aerial status` prints an `ATTENTION:` line under a map that needs a person, with the reason.
- `aerial status --full` prints the whole run record, including every DroneDeploy id.
- A map that says `failed` keeps whatever products it had; nothing is deleted on failure.
- `docs/operations.md` in the module repo has the state table and the failure modes we have met.
