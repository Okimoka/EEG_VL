"""Subject 001 Picard comparison using ordinary MNE-BIDS-Pipeline steps only.

Prepared BIDS sidecars supply geometry, shared bad channels and analysis-event
corrections. ICA uses the native recording-based fixed grid, independently of
those experimental event labels. Stop after classification; review the native
*_proc-ica_components.tsv before any future application step.
"""
from pathlib import Path

_root = Path(__file__).resolve().parent
bids_root = _root / "prepared_native_bids"
deriv_root = _root / "artifacts" / "ica_picard_pilot"
subjects = ["001"]
runs = ["01", "02"]
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
data_type = "eeg"
n_jobs = 2
mne_log_level = "info"  # Native fitting-time and data-loading messages.
# Native caching is enabled by default in this isolated deriv_root.

l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0
raw_resample_sfreq = 250.0
# eeg_reference = "average"  # Default; prepared bads agree across both runs.

task_is_rest = True  # Fixed fitting windows, not a claim of resting-state EEG.
epochs_tmin = 0.0
epochs_tmax = 3.996  # Exactly 1,000 samples at 250 Hz, no shared endpoint.
rest_epochs_duration = 4.0
rest_epochs_overlap = 0.0
baseline = None

spatial_filter = "ica"
ica_reject = "autoreject_local"
# autoreject_n_interpolate = [4, 8, 16]  # Default; no interpolation before ICA.
ica_algorithm = "picard-extended_infomax"
# ica_l_freq = 1.0  # Default.
ica_h_freq = 100.0
# ica_n_components = None  # Default: omit numerically zero PCA dimensions.
# ica_max_iterations = 500  # Default.
ica_use_icalabel = True
ica_use_eog_detection = False
ica_use_ecg_detection = False
random_state = 2026
run_source_estimation = False
decode = False
