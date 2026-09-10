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
pip install "git+ssh://git@github.com/mosscoder/mpg-aerial-pipeline@v0.4.3"
```

The repo is private, so the install goes over SSH with your GitHub key; the `https` form fails.

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
- where the ground control for this flight is, and what it needs to become. DroneDeploy takes ground
  control as a CSV with exactly this header, latitude and longitude in WGS84 degrees, elevation in metres:

  ```
  GCP Label,Latitude,Longitude,Elevation (m)
  SE,46.67531526,-114.00414073,1199.819
  ```

  The Emlid export is a different, much wider table, so it has to be reshaped to those four columns.
  The file then has to be **in the bucket, under `root`**, not on your laptop: the module hands
  DroneDeploy a link to it, so upload the CSV (Cloud Console or `gsutil cp`) and give the config its
  path relative to `root`
- where in DroneDeploy the maps should go: folders and a project, written as a path with ` / ` between
  the levels; the module creates what is missing
- where the products should land in the bucket

## 3. Write the config

Save it as `survey_config.toml` in this folder. One line per decision; the comments say what each is.

```toml
[survey]
id            = "some_project_260101"                 # a unique name for this survey; labels the run records and the log lines
root          = "gs://bucket/surveys/some_project/some_flight"   # bucket prefix holding the raw flight; read only
raw           = "path/under/root/to/DCIM"             # folder under root that holds the DJI flight folders
stem          = "some_project-{yymmdd}"               # name of everything the map produces; {yymmdd} and {date} come from
                                                      #   the map's date, any other {key} from the map block
outputs       = "gs://bucket/surveys/some_project/some_flight/processing/drone_deploy"   # where products and records land
multispectral = true                                  # build the four-band multispectral map from the DJI-calibrated frames
visible       = true                                  # build the RGB orthomosaic and the point cloud

[dronedeploy]
path = "Some Folder / Some Subfolder / Project Name"  # folders / project in DroneDeploy, created if missing;
                                                      #   the plans are named <stem>-visible and <stem>-multispectral

[visible]
gcps = "path/under/root/to/gcps-dronedeploy.csv"      # the GCP CSV IN THE BUCKET, written relative to root (so this one is
                                                      #   gs://bucket/surveys/some_project/some_flight/path/under/root/to/gcps-dronedeploy.csv);
                                                      #   upload it there first. Sent with the visible upload only: DroneDeploy
                                                      #   refuses ground control on multispectral uploads

[maps.flight1]                                        # one block per map; the id is what you use on the command line
                                                      #   (aerial submit --map flight1)
date           = "2026-01-01"                         # the flight date; fills {yymmdd} / {date} above
flight_folders = ["DJI_202601011000_001_mission",     # the DJI flight folders under raw that make up this map, in order;
                  "DJI_202601011030_002_mission"]     #   several when the mission spanned batteries
notes          = "what the next person should know"   # optional
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
