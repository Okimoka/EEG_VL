"""Final cohort ICA: native four-second windows and extended Infomax.

A fresh reproduction includes all 17 subjects, including subject 001.
"""
from pathlib import Path

_root = Path(__file__).resolve().parent
bids_root = _root / "prepared_native_bids"
deriv_root = _root / "artifacts" / "ica_native_cohort"
subjects = [f"{subject:03d}" for subject in range(1, 18)]
runs = ["01", "02"]
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
data_type = "eeg"
n_jobs = 4
mne_log_level = "info"
# Native caching and the loky backend retain their defaults.

l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0
raw_resample_sfreq = 250.0
# eeg_reference = "average"  # Default; prepared bads agree across both runs.

task_is_rest = True
epochs_tmin = 0.0
epochs_tmax = 3.996  # Exactly 1,000 samples at 250 Hz, no shared endpoint.
rest_epochs_duration = 4.0
rest_epochs_overlap = 0.0
baseline = None

spatial_filter = "ica"
ica_reject = "autoreject_local"
# autoreject_n_interpolate = [4, 8, 16]  # Default; no interpolation before ICA.
ica_algorithm = "extended_infomax"
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
