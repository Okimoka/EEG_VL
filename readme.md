# EEG gameplay project

Sana Hafeez, Hadi Ismail, Okan Mazlum · SS 2026

The main report is in `report/report.pdf`. The remaining repository contain the relevant scripts and data concerning the reproduction of this report.


## Files

| Folder/file | Contents |
|---|---|
| `00_prepare_dataset/` | Preparation script and the fixed electrode/coordinate tables |
| `01_bad_channels/` | PyPREP and `review-channels` |
| `configs/` | All mentioned mne-bids-pipeline configs |
| `figures_and_numbers/` | Scripts for the report figures and dataset numbers |
| `figures/` | All standalone images, including unused ones |
| `logs/` | Various process outputs
| `reports/` | Original pipeline HTML reports |
| `report/` | Typst report |
| `analysis/` | Everything relating to the analysis part |
| `requirements.txt`, `setup-environment` | Pinned environment and installation command |

## Setup

Run commands from the submission/repository root:

```bash
./setup-environment
source .venv/bin/activate
```

Download [ds003517 version 1.1.0](https://openneuro.org/datasets/ds003517/versions/1.1.0) and symlink to be visible to scripts

```bash
ln -s /absolute/path/to/ds003517 v1.0.0
```

## Reproduce preprocessing

Prepare the initial dataset. This fixes channel types/counts and copies the supplied corrected BESA coordinates. Original signal files remain unchanged links.

```bash
python 00_prepare_dataset/prepare_bids.py
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

The bundled `electrodes.tsv` and `coordsystem.json` were produced by the earlier preparation script from the dataset's BESA geometry. All 34 recordings share these template positions. The coordinates are expressed in metres, with the head origin and axes defined by the template's nasion and two ear landmarks; they are not individual digitization.

The short `00_prepare_dataset/make_electrodes.py` reproduces the same electrode table using MNE's built-in [BESA reader](https://mne.tools/stable/generated/mne.channels.read_custom_montage.html) and [head-coordinate conversion](https://mne.tools/stable/generated/mne.channels.transform_to_head.html). It downloads the [authors' named BESA template](https://github.com/sccn/dipfit/blob/master/standard_BESA/standard-10-5-cap385.elp), uses its 85 mm radius and preserves the original dataset's 0.01 mm rounding. This is optional: preparation normally just copies the already-generated table.

```bash
python 00_prepare_dataset/make_electrodes.py
```

Reproduce Figure 1 (both subject reports and the 17-subject grand average):

```bash
python .venv/bin/mne_bids_pipeline --config configs/config_first_attempt.py --steps init,preprocessing/data_quality,preprocessing/frequency_filter,preprocessing/make_epochs,preprocessing/ptp_reject,sensor/make_evoked,sensor/group_average
```

Filter at 0.1–100 Hz with a 60 Hz notch, then run all PyPREP checks with seed 2026 and 1,000 RANSAC subsets:

```bash
python .venv/bin/mne_bids_pipeline --config configs/config.py --steps init,preprocessing/data_quality,preprocessing/frequency_filter
python 01_bad_channels/run_pyprep.py
```

PyPREP runs four recordings at a time. It saves flags and review intervals, then prepares `prepared_final_bids/`: each subject gets the union of both runs' flags; ten analysis-event types receive +40 ms; wall/star repeats less than 500 ms after the last retained same-type event receive `EXCLUDED_REPEAT__` labels. Event rows and original values are preserved. Status/start markers are not shifted.

To use the already-saved PyPREP decisions instead of rerunning detection:

```bash
python 00_prepare_dataset/prepare_bids.py --final
```

Inspect channels after recreating the initial filtered data:

```bash
./01_bad_channels/review-channels --subject 001 --run 01
```

Run final full ICA for all 17 subjects:

```bash
python .venv/bin/mne_bids_pipeline --config configs/config_ica_final_full.py --steps preprocessing/data_quality,preprocessing/frequency_filter,preprocessing/fit_ica,preprocessing/find_ica_artifacts
```

This uses 250 Hz, four-second fixed windows, 1–100 Hz ICA training, AutoReject window selection, extended Infomax and ICLabel. It stops before applying exclusions. The recorded run took about two hours. Saved exclusion indices belong to the original models; review any newly fitted models before applying them.

`config_ica_pilot.py` records the earlier event-centred pilot settings. Its reports are included for comparison; reproducing that abandoned custom input runner is outside this minimal package. The controlled Picard test changed the final solver to `picard-extended_infomax`, used subject `001` and a separate output directory. Recorded fitting times: Infomax **528.9 s**, Picard **1006.7 s**.

## Analysis

See [`analysis/README.md`](analysis/README.md) for the short Python scripts, bundled cleaned EEG, saved models and reproduction commands. To regenerate the analysis tables and figure without retraining:

```bash
python analysis/results.py
```

## Numbers and figures

```bash
python figures_and_numbers/dataset_numbers.py
python figures_and_numbers/plot_preprocessing.py
python figures_and_numbers/plot_ica_first10_comparison.py
python figures_and_numbers/plot_eog_examples.py
python figures_and_numbers/plot_repeat_event_example.py
```

The first script reads original events and behavioural logs. The next two use the bundled small outputs. EOG examples need the original dataset; the repeated-wall example also needs `prepared_final_bids/`. Images go to `figures/`; numeric summaries go to `logs/`.

Compile with the external Typst CLI:

```bash
typst compile --root . report/report.typ report/report.pdf
```

## Working copy and upload

The submission scripts are now independent copies, so simplifying them does not modify the working project's earlier runs. While editing the main project report, refresh this copy and its images with:

```bash
python submission/update-report.py
```

That command runs from the parent working project. It only updates the submission copy's image paths and run names. Standalone HTML reports and small outputs may still be linked to their originals. Before uploading, copy this folder to a **new** destination to replace those links with ordinary files:

```bash
cp -RL submission github_submission
```

Do this before adding a local dataset or environment inside `submission/`, so those are not copied into the upload folder.
