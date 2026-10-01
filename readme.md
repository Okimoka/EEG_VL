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

Prepare the initial dataset (fix channel and electrodes sidecars)

```bash
python 00_prepare_dataset/prepare_bids.py
```

Filter at 0.1–100 Hz with a 60 Hz notch, then run all PyPREP checks with seed 2026 and 1,000 RANSAC subsets:

```bash
python .venv/bin/mne_bids_pipeline --config configs/config.py --steps init,preprocessing/data_quality,preprocessing/frequency_filter
python 01_bad_channels/run_pyprep.py
```

This should output `prepared_final_bids/`: each subject gets the union of both runs' flags; ten analysis-event types receive +40 ms; wall/star repeats less than 500 ms after the last retained same-type event receive `EXCLUDED_REPEAT__`

Inspect channels after recreating the initial filtered data:

```bash
./01_bad_channels/review-channels --subject 001 --run 01
```

Run final full ICA for all 17 subjects:

```bash
python .venv/bin/mne_bids_pipeline --config configs/config_ica_final_full.py --steps preprocessing/data_quality,preprocessing/frequency_filter,preprocessing/fit_ica,preprocessing/find_ica_artifacts
```

This uses 250 Hz, four-second fixed windows, 1–100 Hz ICA training, AutoReject window selection, extended Infomax and ICLabel. It stops before applying exclusions. The recorded run took about two hours.

## Analysis

See [`analysis/README.md`](analysis/README.md) for the short Python scripts, bundled cleaned EEG, saved models and reproduction commands. To regenerate the analysis tables and figure without retraining:

```bash
python analysis/results.py
```