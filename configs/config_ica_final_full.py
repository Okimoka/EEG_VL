from pathlib import Path

_root = Path(__file__).resolve().parents[1]
bids_root = _root / "prepared_final_bids"
deriv_root = _root / "logs" / "ica_final_full"
subjects = [f"{subject:03d}" for subject in range(1, 18)]
runs = ["01", "02"]
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
data_type = "eeg"
n_jobs = 4
mne_log_level = "info"

l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0
raw_resample_sfreq = 250.0

task_is_rest = True
epochs_tmin = 0.0
epochs_tmax = 3.996  # Exactly 1,000 samples at 250 Hz
rest_epochs_duration = 4.0
rest_epochs_overlap = 0.0
baseline = None

spatial_filter = "ica"
ica_reject = "autoreject_local"
ica_algorithm = "extended_infomax"
ica_h_freq = 100.0
ica_use_icalabel = True
ica_use_eog_detection = False
ica_use_ecg_detection = False
random_state = 2026
run_source_estimation = False
decode = False
