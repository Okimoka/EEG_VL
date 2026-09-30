"""Apply reviewed ICA, then clean event-related epochs with the native pipeline.

Never run fit_ica or find_ica_artifacts with this analysis configuration.
EEG_POSTICA_AVERAGING_GROUP only handles missing gambling conditions when
creating native evoked reports; preprocessing always uses the default 'all'.
"""
from pathlib import Path
import os

_root = Path(__file__).resolve().parent
bids_root = _root / "prepared_native_bids"
deriv_root = _root / "artifacts" / "postica"
subjects = [f"{number:03d}" for number in range(1, 18)]
runs = ["01", "02"]
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
data_type = "eeg"
n_jobs = 4
mne_log_level = "info"

l_freq = 0.1
h_freq = 40.0
notch_freq = 60.0
# raw_resample_sfreq = None  # Default: preserve the original 500 Hz EEG/EOG.
# eeg_reference = "average"  # Default: use the same good-channel set as ICA.
# task_is_rest = False  # Default: actual prepared events now define epochs.
conditions = ["ODDBALL STANDARD", "ODDBALL RARE", "GAMBLING WIN", "GAMBLING LOSS",
              "PLAYER_CRASH_WALL", "COLLECT_STAR", "MISSILE_HIT_ENEMY",
              "PLAYER_CRASH_ENEMY", "COLLECT_AMMO", "SHOOT_BUTTON"]
epochs_tmin = -0.2
epochs_tmax = 0.8
baseline = (-0.2, 0.0)

spatial_filter = "ica"
ica_use_icalabel = True  # Apply the same average-reference convention.
ica_use_eog_detection = False
ica_use_ecg_detection = False
ica_algorithm = "extended_infomax"
ica_h_freq = 100.0  # Describes the saved fit; there is no fitting in this run.
ica_reject = "autoreject_local"
reject = "autoreject_local"
# autoreject_n_interpolate = [4, 8, 16]  # Default; actual repairs occur after ICA.
random_state = 2026
contrasts = [("ODDBALL RARE", "ODDBALL STANDARD"), ("GAMBLING WIN", "GAMBLING LOSS")]
run_source_estimation = False
decode = False

_group = os.environ.get("EEG_POSTICA_AVERAGING_GROUP", "all")
_missing_gambling = {"002", "003", "005"}
if _group == "gambling_present":
    subjects = [subject for subject in subjects if subject not in _missing_gambling]
elif _group == "gambling_absent":
    subjects = sorted(_missing_gambling)
    conditions = [condition for condition in conditions if not condition.startswith("GAMBLING")]
    contrasts = contrasts[:1]
elif _group != "all":
    raise ValueError(f"Unknown averaging group: {_group}")
