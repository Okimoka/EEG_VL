# EOG treatment

VEOG/HEOG samples are preserved. They do not enter PyPREP, the EEG average reference, ICA fitting or automatic component selection. ICLabel and visual EEG-component review determine the ICA exclusions.

`reference/qc/eog_metrics.tsv` records full-recording offsets, excursions and peak-to-peak values under the loader's EEGLAB scaling convention. `scripts/audit_dataset.py` regenerates these values. `reference/eog_channel_audit.tsv` and `scripts/audit_eog.py` use the middle 120 s of each recording, a 0.5–15 Hz Butterworth filter and five-second edge removal to compare each auxiliary channel with Fp1−Pz. This diagnostic is separate from pipeline filtering. Subjects 003, 004 and 006 had weak VEOG/frontal association and small filtered variation in the inspected intervals.

The large, heterogeneous auxiliary values do not establish a universal physiological gain correction. Native `.fdt` values and MNE values agree under the expected 10⁻⁶ conversion. Subject 014's suspected exchanged EOG roles remain a visual suspicion, not a verified wiring diagnosis. No relabelling, rescaling or interpolation of EOG is applied.
