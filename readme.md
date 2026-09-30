# EEG gameplay reproduction — preprocessing submission

Sana Hafeez, Hadi Ismail, Okan Mazlum · SS 2026

Preprocessing for Cavanagh and Castellanos (2016), [*Identification of canonical neural events during continuous gameplay of an 8-bit style video game*](https://doi.org/10.1016/j.neuroimage.2016.02.075), using [OpenNeuro ds003517](https://openneuro.org/datasets/ds003517/versions/1.1.0).

This folder is a standalone repository: upload **its contents**, rather than the surrounding working project. It contains the selected final preprocessing workflow, accepted decisions and compact evidence for the report. Classification/transfer analysis will be added separately. Original EEG, large derivatives, HTML reports, caches, course materials and superseded pilot pipelines are not bundled. They are either downloaded or regenerated locally.

## Contents

| Path | Purpose |
|---|---|
| `config.py` | Preliminary 0.1–100 Hz filtering at 500 Hz for channel detection. |
| `config_ica_native_cohort.py` | Final four-second-window ICA for **all 17 subjects**, including 001. |
| `config_postica.py` | 0.1–40 Hz event epochs, reviewed ICA application, fresh AutoReject and baseline. |
| `config_broadband.py` | Optional 0.1–100 Hz continuous ICA-cleaned output at 500 Hz. |
| `scripts/` | Metadata preparation, standalone PyPREP, small pipeline launch/staging helpers, validation and report figures. |
| `01_review_bad_channels/` | The channel inspection browser used during review. |
| `reference/` | Accepted detector JSONs, 17 small ICA models and their reviewed component TSVs, summary tables, a compressed trial ledger and compact plot inputs. These are recorded results, not a workspace for new runs. |
| `report/` | Standalone snapshot of the introduction/preprocessing chapter, its figures and bibliography. Existing inline notes/placeholders are retained; this is not the complete analysis report. |
| `docs/` | Short geometry, EOG and provenance explanations. |
| `checks/sub007/` | Optional, isolated reproduction of the reported rejection-policy check. |
| `requirements.txt` | Recorded Python package versions and the exact course-pipeline commit. |

The workflow is:

1. Prepare a metadata-corrected BIDS view; original samples remain unchanged.
2. Filter at 0.1–100 Hz with a 60 Hz notch; run PyPREP separately (seed 2026, 1,000 RANSAC subsets).
3. Prepare the subject-wise union of bad channels, +40 ms event correction and retained-first 500 ms wall/star rule.
4. Fit one extended Infomax ICA per subject across both runs, at 250 Hz, using non-overlapping four-second windows and fresh local AutoReject selection. No interpolation or baseline correction is applied to ICA training data.
5. Review ICLabel proposals and apply accepted ICA exclusions to newly filtered 500 Hz, 0.1–40 Hz event epochs. A second AutoReject fit repairs/rejects trials, followed by a −0.2…0 s baseline.
6. Validate trial counts and create matched before/after cleaning ERPs at Pz (oddball) and Cz (gambling).

## Environment and dataset

Use Linux, Python **3.12** (recorded: 3.12.3), Git and a Python virtual environment. The Qt channel browser needs a graphical desktop. Typst **0.15.1** is optional for compiling the report. The commands below run from this repository's root.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/setup_pipeline.py
python scripts/setup_pipeline.py --check
python -m unittest discover -s scripts/tests -v
```

The setup helper applies/verifies the recorded upstream compatibility change to the pinned course fork. It is needed when recreating the environment; it changes no data. Do not silently substitute a newer pipeline or different package versions when comparing numerical results.

Download the complete dataset, including signal payloads (`*.fdt`), sidecars and `code/Logs/`, from OpenNeuro. Put it at `v1.0.0/`, or create a local link:

```bash
ln -s /absolute/path/to/downloaded/ds003517 v1.0.0
python scripts/verify_inputs.py --data v1.0.0
```

The directory name `v1.0.0` is historical; the supplied `dataset_description.json` identifies **dataset DOI version 1.1.0**. The source manifest checks the actual files used here, including EEG bytes, rather than relying on the folder name. A Git/DataLad checkout with missing annex payloads is not enough. In the original working project, the corresponding link is `ln -s ../v1.0.0 v1.0.0` from inside `submission/`.

The original data occupy about 6.5 GB; regenerated derivatives and reports need considerably more space. Allow roughly 80 GB of free disk space for a full reproduction with both output branches. The main configurations use four workers; reduce `n_jobs` and PyPREP's `--jobs` if memory is limited. The launcher limits each worker's numerical-library threads to avoid oversubscription.

## Full reproduction from raw data

### 1. Initial preparation, frequency filtering and PyPREP

```bash
python scripts/prepare_bids.py --source v1.0.0 --destination prepared_bids
python scripts/run_pipeline.py --log artifacts/mne-bids-pipeline/native_pipeline.log --config config.py --steps init,preprocessing/data_quality,preprocessing/frequency_filter
python scripts/run_pyprep.py --jobs 4
python scripts/channel_policy.py
python scripts/prepare_bids.py --source v1.0.0 --destination prepared_native_bids --analysis-ready
```

Preparation corrects EEG/EOG types and counts, units and electrode coordinates. The second invocation also writes the common bad-channel set and analysis-event metadata. It always reads the **original** dataset. It keeps all event rows and preserves original onsets/labels in extra columns; excluded repeats receive an `EXCLUDED_REPEAT__` prefix. Never add another 40 ms downstream. EOG samples and the short subject-008 gameplay tail are retained.

The adopted policy accepts every PyPREP threshold flag without overrides. `channel_policy.py` binds the exported policy to the detector files and reports differences from the accepted reference. Changes in channel selection require renewed review before using the saved reference ICA models.

For optional visual channel review:

```bash
python 01_review_bad_channels/review_store.py
./01_review_bad_channels/review-channels --subject 001 --run 01
```

See [the browser instructions](01_review_bad_channels/readme.md). The workbook is a local inspection aid; the reproduction policy remains **accept all PyPREP flags**. Editing its overrides does not silently change the exported policy.

### 2. Fit and review ICA

```bash
python scripts/run_pipeline.py --log artifacts/ica_native_cohort/native_pipeline.log --config config_ica_native_cohort.py --steps preprocessing/data_quality,preprocessing/frequency_filter,preprocessing/fit_ica,preprocessing/find_ica_artifacts
python scripts/summarize_ica.py
```

Open `artifacts/ica_native_cohort/sub-XXX/eeg/sub-XXX_task-ContinuousVideoGamePlay_report.html`. These are the normal pipeline reports with ICLabel output and ICA decompositions. Review maps, spectra and time courses, then update the native `*_proc-ica_components.tsv` files: `component` is a **zero-based** index; `status` is `bad` or `good`; document changes in `status_description`.

The original run had 178 automatic proposals and 15 additional manual exclusions. The accepted tables and reasons are in `reference/ica/`, including `manual_additions.tsv`. **Do not copy those component numbers onto a newly fitted model without matching/reviewing it.** ICA order and separation can change with numerical libraries or inputs. The shorter replay route below uses the actual reviewed models and requires no new component decisions.

After reviewing the fresh fit:

```bash
python scripts/prepare_ica_application.py --models refit --reviewed
```

This stages models and decisions in a separate derivative directory and builds the eligible-event ledger. It does not refit ICA or modify EEG. Do not rerun `find_ica_artifacts` after editing the component TSVs, because it rewrites automatic proposals.

### 3. Apply ICA and clean the analysis epochs

```bash
python scripts/run_pipeline.py --log artifacts/postica/native_pipeline.log --config config_postica.py --steps preprocessing/frequency_filter,preprocessing/make_epochs,preprocessing/apply_ica,preprocessing/ptp_reject
EEG_POSTICA_AVERAGING_GROUP=gambling_present python scripts/run_pipeline.py --log artifacts/postica/native_pipeline.log --config config_postica.py --steps sensor/make_evoked
EEG_POSTICA_AVERAGING_GROUP=gambling_absent python scripts/run_pipeline.py --log artifacts/postica/native_pipeline.log --config config_postica.py --steps sensor/make_evoked
python scripts/analyze_postica.py
```

The two averaging calls handle subjects 002, 003 and 005, which have no gambling outcomes; all subjects receive the same preprocessing. `analyze_postica.py` checks the reference and ICA application, unchanged EOG, corrected event identities, final trial selection and measured repairs. It creates counts, the trial ledger and subject averages on identical final retained trials. It does not perform new rejection or fitting.

Key outputs in `artifacts/postica/sub-XXX/eeg/`:

- `*_proc-clean_epo.fif`: final −0.2…+0.8 s, 500 Hz, 0.1–40 Hz epochs with reviewed ICA removal, local AutoReject repairs/rejection and baseline. Global bad channels remain marked.
- `*_run-XX_proc-clean_raw.fif`: continuous ICA-cleaned 0.1–40 Hz data. **No trial-specific AutoReject repairs or baseline subtraction.**
- `*_report.html`: native preprocessing/application reports.

Subject averages in `artifacts/postica/evokeds/` additionally interpolate globally bad channels for display. The ERP figure's lower row includes both ICA removal and final AutoReject repairs; the separate ICA-only comparison omits those repairs from both curves.

### Shorter replay using the accepted decisions

After environment setup and input verification, the saved detector decisions and reviewed ICA models can reproduce the final processing without refitting them:

```bash
python scripts/prepare_bids.py --source v1.0.0 --destination prepared_bids
python scripts/channel_policy.py --reference
python scripts/prepare_bids.py --source v1.0.0 --destination prepared_native_bids --analysis-ready
python scripts/prepare_ica_application.py --models reference
```

Then run **step 3** above. Model and decision-file hashes, good-channel order and geometry are checked before application. This route regenerates application reports; recreating the original training reports and ICLabel tables requires **step 2**. The original model files retain their automatic exclusion list internally; the paired component-status TSVs supply the complete **193 reviewed exclusions**, as in the native pipeline.

### Optional broadband continuous output

Use the same reviewed models, with either `--models reference` or `--models refit --reviewed` consistently with the chosen route:

```bash
python scripts/prepare_ica_application.py --branch ica_broadband --models reference
python scripts/run_pipeline.py --log artifacts/ica_broadband/native_pipeline.log --config config_broadband.py --steps preprocessing/frequency_filter
python scripts/apply_broadband_ica.py
```

This saves `artifacts/ica_broadband/sub-XXX/eeg/*_proc-clean_raw.fif`: **0.1–100 Hz, 500 Hz**, 60 Hz notch, common average reference and the same ICA exclusions. The small adapter calls the pipeline's native continuous-ICA function because this branch needs no event epochs. It does not fit a second ICA. Any later long-epoch/time–frequency analysis must define its own epoch lengths, margins and cleaning policy; repairs from the short ERP epochs cannot be transferred automatically.

## Report figures and checks

To redraw the main report plots from the compact accepted reference data, without raw EEG or expensive processing:

```bash
python scripts/verify_inputs.py
python scripts/make_report_preprocessing_figures.py --reference
python scripts/make_report_erp_stages.py --reference
```

The same scripts without `--reference` use recomputed `artifacts/` results. The first then requires the preliminary filtered recording from step 1 and `summarize_ica.py`; the second requires `analyze_postica.py`. In the shorter replay route, use `--reference` for the recorded ICA-training/channel plots.

The compact ERP files contain only Pz/Cz for the four exemplar conditions, after the original display interpolation/baseline; they are plot inputs, **not classifier-ready epochs**. Blink snippets contain the exact displayed four-second before/after traces. The original initial-inspection screenshot and example component-property image are retained as static figures. The three channel-example placeholders and any outstanding inline report notes are preserved from the supplied report, not fabricated results.

Additional read-only audits after preparation:

```bash
python scripts/audit_dataset.py
python scripts/audit_eog.py
python scripts/audit_requested_numbers.py
python scripts/audit_game_alignment.py
```

The latter two use behavioural logs and ICA component TSVs; run them after fitting/reviewing ICA. They reproduce task durations, recording exceptions, eye-component counts and the sub-003/sub-008 log alignment. For reference replay, the original audit tables are already in `reference/report_notes_revision/`.

Optional checks discussed in the report:

- `python scripts/diagnose_native_ica_rejection.py` after the full fit replays AutoReject for 005/014 and compares its selection with a 500 µV rule. It performs no new ICA fit. Its task annotations come from the archived post-hoc window audit and do not affect selection.
- The Picard solver comparison is described in [docs/reproduction_notes.md](docs/reproduction_notes.md). Earlier event-window/custom-block pilots remain historical checks. The small subject-007 sensitivity scripts are separate under `checks/sub007/`; their outputs never replace main results (see the same notes).

Compile the preprocessing snapshot with:

```bash
typst compile --root . report/preprocessing.typ report/preprocessing.pdf
```

## Recorded reference checks

| Check | Accepted result |
|---|---:|
| Subjects / recordings | 17 / 34 |
| PyPREP channel–recording flags | 68 |
| Ventral flags / other-electrode flags | 21/136 (15.4%) / 47/2,006 (2.3%) |
| ICA fitting windows retained | 9,714 / 11,655 |
| ICA components / automatic exclusions / reviewed exclusions | 1,003 / 178 / 193 |
| Eligible analysis-event anchors | 36,258 |
| Available analysis epochs | 36,244 (13 boundary windows and one out-of-recording event omitted) |
| Final epochs / rejected epochs | 32,195 / 4,049 |
| Retained epochs with measured channel repair | 28,316 (88.0%) |

These are recorded results, not tuning targets. Refit differences must be investigated and reported. High trial loss in 004/007 and frequent repairs remain limitations; the package does not silently replace those subjects or relax their thresholds.

The original cohort fit took roughly two hours on the project machine; post-ICA processing took about 28 minutes. PyPREP with 1,000 RANSAC subsets is a separate substantial job. Hardware, worker count and caching affect runtimes. Stages save logs and use native pipeline caching; after a failure, resolve the error and repeat that stage rather than deleting previous outputs.

`.gitignore` keeps original data, prepared views, environments and generated `artifacts/` out of Git. Keep `reference/`: it is intentionally small and binds the reported decisions to the actual models. Dataset authorship and CC0 terms remain in the original OpenNeuro metadata; this repository does not relicense upstream software or redistribute the paper/course materials.
