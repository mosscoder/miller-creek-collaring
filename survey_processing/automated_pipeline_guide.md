# Processing the 2026-09-08 collaring flight with the automated pipeline

The `aerial` module takes a flight's raw frames from the survey bucket, builds the DroneDeploy uploads
in the cloud (the multispectral frames get DJI's radiometric calibration on the way), creates the maps in
DroneDeploy, waits for them to process, lands the exports back in the bucket, and then registers the
multispectral map to the visible map on a GPU. Nothing passes through your laptop except one config file
and a few commands. The module lives at https://github.com/mosscoder/mpg-aerial-pipeline; its
`docs/survey_config.md` explains every config key and `docs/operations.md` every moving part.

This is a homework assignment: the config below is a template, and you fill it in from what you know
about the flight. The preflight command tells you what is wrong until it is right.

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

## 2. Find out what the config needs to say

The flight was uploaded under `gs://mpg-aerial-survey/surveys/miller_collaring/` in a folder named by
its date. Before writing anything, work out:

- which DJI flight folders make up the mission, and which folders under `DCIM` are not part of it
  (test captures, panel shots). `aerial discover --raw gs://…/DCIM` drafts map blocks from a bucket
  listing and counts the frames in each folder, which helps here
- where the ground control for this flight is, and what it needs to become: DroneDeploy takes ground
  control as a CSV with the header `Label,Latitude,Longitude,Elevation` in WGS84, and the module reads
  it from a path inside the bucket that the config names
- where in DroneDeploy the maps should go: folders and a project, written as a path with ` / ` between
  the levels; the module creates what is missing
- where the products should land in the bucket

## 3. Write the config

Save it as `survey_config.toml` in this folder. One line per decision; the comments say what each is.

```toml
[survey]
id            = "..."                    # a unique name for this survey; labels the run records and the log lines
root          = "gs://..."               # bucket prefix holding the raw flight; read only
raw           = "..."                    # folder under root that holds the DJI flight folders (e.g. "data_collection/DCIM")
stem          = "..."                    # name of everything the map produces; may use placeholders from the map block,
                                         #   e.g. "miller_collaring-{yymmdd}" when the map has a date
outputs       = "gs://..."               # where products and records land (and the working copy of this file)
multispectral = true                     # build the four-band multispectral map from the DJI-calibrated frames
visible       = true                     # build the RGB orthomosaic and the point cloud

[dronedeploy]
path = "..."                             # folders / project in DroneDeploy, e.g. "Side Projects / Miller Creek";
                                         #   the plans are named <stem>-visible and <stem>-multispectral

[visible]
gcps = "..."                             # the GCP CSV, relative to root; sent with the visible upload only
                                         #   (DroneDeploy refuses ground control on multispectral uploads)

[maps.<id>]                              # one block per map; the id is what you use on the command line (aerial submit --map <id>)
date           = "YYYY-MM-DD"            # the flight date; fills {yymmdd} / {date} in the templates above
flight_folders = ["...", "..."]          # the DJI flight folders under raw that make up this map, in order
notes          = "..."                   # anything the next person should know (optional)
```

Everything else (frame patterns, export layers, projection, the calibration rule) is the module's default.

## 4. Run it

```bash
aerial submit --dry-run     # preflight: counts the frames, resolves the DroneDeploy path, checks the GCP file; creates nothing
aerial submit               # builds both uploads in the cloud, creates the two plans, posts the transfers
aerial status               # one line per map and data type, from the run records in the bucket
```

Iterate on the dry run until it is clean. `submit` then adds this config to the watch list of the poll
job, which runs every 15 minutes in Cloud Run and drives everything from there: it waits for DroneDeploy
to fetch the uploads, follows the two maps through processing, requests the exports, lands them in the
bucket, and finally runs the registration on the GPU job. You can close the laptop after `submit`.

One step needs a person. Ground control makes DroneDeploy tag the targets in the images and then wait
for someone to review the tags: open the visible plan in the DroneDeploy app, check the tags, and press
Continue to Processing. `aerial status` says `gcp=pending` while it waits and flags the map after two
hours. The multispectral map has no ground control and processes on its own.

## 5. Read the result

Run `aerial status` whenever you like. A map goes `submitted` → `transferred` → `processing` →
`processed` → `exporting` → `exported` → `done`; the multispectral line shows `steps: register` at the
end. Under `outputs` you will find:

```
<stem>-visible.tif               RGBA orthomosaic, GCP-rectified, EPSG 6514
<stem>-pointcloud.las
<stem>-multispectral.tif         Red, Green, NIR, RedEdge + alpha, registered to the visible map
<stem>-*.manifest.json           every frame uploaded, md5s, the calibration terms
<stem>-*.run.json                the state of each map, step by step
registration_results/<stem>/     quality.tif, qa.json, quicklook.png, report.html
survey_config.toml               the working copy of your file
```

`report.html` is self-contained: download it and open it in a browser. `quality.tif` is one small raster
on a 64-pixel cell grid with the measured shifts, match coverage and the agreement with the visible map
before and after; `qa.json` has the same numbers plus timing.

Expect about a day end to end, most of it DroneDeploy's queue, and on the order of ten dollars of
cost, almost all of it DroneDeploy pulling the uploads out of the bucket.

## If something looks wrong

- `aerial status` prints an `ATTENTION:` line under a map that needs a person, with the reason.
- `aerial status --full` prints the whole run record, including every DroneDeploy id.
- A map that says `failed` keeps whatever products it had; nothing is deleted on failure.
- `docs/operations.md` in the module repo has the state table and the failure modes we have met.
