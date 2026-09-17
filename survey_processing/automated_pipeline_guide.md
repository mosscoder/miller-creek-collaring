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

**Your own key.** Everything the pipeline does in the cloud runs as one identity, a service account,
and the laptop side runs as that same account through a key file. Ask for your own key on it rather
than copying someone else's: same identity and permissions, but yours to revoke, and the audit log can
tell whose laptop a run came from. Keep the file outside any repo; its path goes in the `.env` below.

Python 3.10 or newer, in whatever environment you like:

```bash
pip install "git+ssh://git@github.com/mosscoder/mpg-aerial-pipeline"
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
- where the ground control file will live in the bucket; section 4 covers the file itself
- where in DroneDeploy the maps should go: folders and a project, written as a path with ` / ` between
  the levels; the module creates what is missing
- where the products should land in the bucket

## 3. Write the config

Save it as `survey_config.toml` in this folder. One line per decision; the comments say what each is.

```toml
[survey]
project       = "some_project"                        # a free key; fills {project} wherever it appears below
date          = "260101"                              # the flight date, in the form you want to see in names; fills {date}
agl           = 50                                    # planned height above ground, metres (the mission's setting: 25, 50 or 125);
                                                      #   the export resolution comes from it and a coarser export is refused
root          = "gs://bucket/surveys/{project}/{date}"        # bucket prefix holding the raw flight; read only
raw           = "path/under/root/to/DCIM"             # folder under root that holds the DJI flight folders
outputs       = "{root}/processing/drone_deploy"      # where products and records land
anchor        = "gs://bucket/surveys/{project}/earlier_flight/processing/drone_deploy/earlier-visible.tif"
                                                      # a tif in the bucket: BOTH this flight's visible and multispectral maps
                                                      #   are registered to it once they land. Leave the line out and the
                                                      #   multispectral map registers to this flight's own visible map instead
multispectral = true                                  # build the four-band multispectral map from the DJI-calibrated frames
visible       = true                                  # build the RGB orthomosaic and the point cloud

[dronedeploy]
path   = "Some Folder / Some Subfolder / Project Name"  # folders / project in DroneDeploy, created if missing; the plans are
                                                        #   named <project>_<date>-visible and <project>_<date>-multispectral
emails = ["you@mpgranch.com"]                           # optional: DroneDeploy emails these the link when an export finishes

[visible]
gcps = "path/under/root/to/gcps-dronedeploy.csv"      # the GCP CSV IN THE BUCKET, written relative to root (so this one is
                                                      #   gs://bucket/surveys/some_project/260101/path/under/root/to/gcps-dronedeploy.csv);
                                                      #   upload it there first. Sent with the visible upload only: DroneDeploy
                                                      #   refuses ground control on multispectral uploads

[maps.flight1]                                        # one block per map; the id is what you use on the command line
                                                      #   (aerial submit --map flight1)
flight_folders = ["DJI_202601011000_001_mission",     # the DJI flight folders under raw that make up this map, in order;
                  "DJI_202601011030_002_mission"]     #   several when the mission spanned batteries
notes          = "what the next person should know"   # optional
```

Everything the map produces is named `<project>_<date>`; write an `id` line under `[survey]` to name it otherwise.
Everything else (the aircraft profile, frame patterns, export layers, projection, the calibration rule) is the
module's default.

## 4. Ground control

We treat every point on this flight as a checkpoint: the frames carry RTK positions, so the points measure
the map rather than bend it. DroneDeploy reads the CSV below (WGS84 degrees, metres) and takes any label
containing the word `checkpoint`, lowercase, as a checkpoint; a label without it would act as a ground
control point. Upload the file into the bucket under `root` and name it in `[visible] gcps`.

```
GCP Label,Latitude,Longitude,Elevation (m)
1_checkpoint,46.67531526,-114.00414073,1199.819
2_checkpoint,46.67612044,-114.00302918,1201.455
```

## 5. Run it

```bash
aerial submit --dry-run     # preflight: counts the frames, resolves the DroneDeploy path, reads the GCP file and shows how
                            #   DroneDeploy will treat each point (a bad header or row is refused); creates nothing
aerial submit               # builds both uploads in the cloud, creates the two plans, posts the transfers
aerial status               # one line per map and data type, from the run records in the bucket
aerial board --open         # the survey board: the same records as one page in the browser, rewritten every 5 minutes
aerial refetch --map <id>   # fetch a landed product again under the current config: downloads the matching DroneDeploy export, renders only if none; re-registers
```

Iterate on the dry run until it is clean. `submit` then adds this config to the watch list of the poll
job, which runs every 5 minutes in Cloud Run and drives everything from there: it waits for DroneDeploy
to fetch the uploads, follows the two maps through processing, requests the exports, lands them in the
bucket, and finally runs the registration on the GPU job. You can close the laptop after `submit`.

One step needs a person. When the visible plan reaches its tag review in the DroneDeploy app, check the
tags, confirm every point shows as a checkpoint, and press Continue to Processing. `aerial status` says
`gcp=pending` while it waits and flags the map after two hours.

## 6. Read the result

Run `aerial status` whenever you like. A map goes `submitted` → `transferred` → `processing` →
`processed` → `exporting` → `exported` → `done`; with an anchor, both lines show `steps: register` at the
end. Under `outputs` you will find:

```
<name>-visible.tif               RGBA orthomosaic at the camera's GSD for the flight height, EPSG 6514
<name>-pointcloud.las
<name>-multispectral.tif         Red, Green, NIR, RedEdge + alpha, registered to the anchor
<name>-*.manifest.json           every frame uploaded, md5s, the calibration terms
<name>-*.run.json                the state of each map, step by step, and the export resolution asked for and delivered
registration_results/<name>-<type>/   quality.tif, qa.json, quicklook.png, report.html
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
- `aerial board --open` opens the survey board: a banner with what needs a person, a progress bar, the deliverables, and one timeline per map with the time of every step. The poller rewrites it every 5 minutes; `--refresh` rewrites it now.
- A map that says `failed` keeps whatever products it had; nothing is deleted on failure.
- `docs/operations.md` in the module repo has the state table and the failure modes we have met.
