# Reproduction boundaries and optional solver check

The submitted code keeps the scientific settings of the selected final workflow. Packaging changes remove dependencies on discarded pilot directories: the ICA config now includes subject 001, channel-policy export directly implements the adopted all-flags rule, and model staging checks the supplied small models instead of requiring historical HTML reports and validation transactions. The pipeline launcher only adds logging and thread limits. Filtering, resampling, referencing, ICA, event epoching, AutoReject and native reports remain MNE-BIDS-Pipeline operations.

Preparation is custom: corrected channel metadata/geometry, shared bad-channel sets, +40 ms event timing and the 500 ms eligibility labels. PyPREP and its inspection browser are separate. Postprocessing scripts validate native outputs and build the report's matched ERPs and diagnostics. Continuous-only broadband application calls the pinned pipeline's existing function through a small adapter.

`reference/` stores original selected results, with hashes in `file_manifest.json`. The electrode Info file, blink snippets and compressed Pz/Cz evokeds are compact extracts for figures; they do not replace full analysis epochs. `trial_ledger.tsv.gz` contains original event-row identifiers, inclusion and measured repair counts. `window_diagnostics.tsv` retains the post-hoc task annotations used to contextualize ICA retention; those inferred task labels were never used to choose the native fixed grid.

For a new environment, input verification, pinned dependencies, fixed seed and thread limits improve reproducibility but do not guarantee bitwise-identical floating-point optimizers on every platform. Fresh ICA indices require review. The accepted reference models and paired TSVs provide an unambiguous replay of the reviewed component decisions.

## Picard comparison (optional)

After reproducing regular Infomax for subject 001 in the cohort directory:

```bash
python scripts/run_pipeline.py --log artifacts/ica_picard_pilot/native_pipeline.log --config config_ica_picard_pilot.py --steps preprocessing/data_quality,preprocessing/frequency_filter,preprocessing/fit_ica,preprocessing/find_ica_artifacts
python scripts/compare_solvers.py
```

The comparison checks that retained training samples match exactly, then finds one-to-one component-map matches and compares automatic exclusions. Native logs contain fit runtimes. The recorded comparison was 16 min 47 s for Picard and 8 min 49 s for regular Infomax, with 57/60 matched maps above 0.998 correlation and identical matched automatic exclusions. `reference/checks/picard_matched_components.tsv` preserves the original matches. The one-subject result is not a general speed benchmark.

Earlier event-centred and custom task-block pilots explain the methodological history in the working report, but are not dependencies of the selected final workflow and are not packaged as additional preprocessing stages. Their entire large archive is excluded. The 005/014 counterfactual rejection summary and 007 policy-comparison counts are retained under `reference/checks/` as supporting historical results; they did not replace the selected fits or final trial cleaning. Reproducing every abandoned pilot would require that separate working archive. The specific subject-007 rejection experiment mentioned in the preprocessing chapter can be rerun after `analyze_postica.py` with:

```bash
python checks/sub007/investigate.py
python checks/sub007/compare.py
python checks/sub007/spatial_checks.py
```

It first verifies the native AutoReject selection, then compares prespecified 75/100 µV screens on the same post-ICA data. Expect roughly 30–40 minutes, depending on hardware. All exploratory files go to `artifacts/checks/sub007/`; no main ICA, trial selection or report result is overwritten. The inspection scripts produce evidence; they do not automatically choose a policy based on ERP/classifier performance.

## Packaging verification

See `reference/package_validation.json` for checks actually executed while preparing this submission. These include metadata unit tests, a clean-directory preparation/configuration smoke check, native ICA application on a short recording copy, regeneration of reference plots and Typst compilation. They are not a claim that the full PyPREP/ICA/AutoReject cohort was rerun during packaging or that a fresh network installation was tested.
