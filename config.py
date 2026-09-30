"""Initial, pre-review configuration for the installed course pipeline fork.

Run only init, data_quality and frequency_filter until channels are reviewed.
ICA settings describe the planned first experiment; no ICA result is implied.
"""
from pathlib import Path

_root = Path(__file__).resolve().parent
bids_root = _root / "prepared_bids"
deriv_root = _root / "artifacts" / "mne-bids-pipeline"
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
data_type = "eeg"
n_jobs = 3

# subjects = "all"
# runs = "all"  # 01: OddGamble; 02: Axon. Same frequency filtering for both.
# eeg_template_montage = None  # Use supplied electrodes.tsv through MNE-BIDS.
# eeg_reference = "average"  # Applied before ICLabel/ICA; frequency-only FIFs retain CPz.
# drop_channels = []  # Inspect the four ventral electrodes, do not discard a priori.
# raw_resample_sfreq = None  # Keep original 500 Hz and event precision.
l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0

conditions = ["ODDBALL STANDARD", "ODDBALL RARE", "GAMBLING WIN", "GAMBLING LOSS",
              "PLAYER_CRASH_WALL", "COLLECT_STAR", "MISSILE_HIT_ENEMY", "PLAYER_CRASH_ENEMY",
              "COLLECT_AMMO", "SHOOT_BUTTON"]
epochs_tmin = -2.0
epochs_tmax = 2.0
baseline = (-0.2, 0.0)
spatial_filter = "ica"
ica_reject = "autoreject_local"
ica_algorithm = "extended_infomax"
# ica_l_freq = 1.0
ica_h_freq = 100.0
ica_use_icalabel = True
# ica_n_components = None  # Rank-aware PCA; do not force 63 components.
ica_use_eog_detection = False  # Unresolved EOG gain; review each channel before enabling.
# VEOG/HEOG remain available for separate, normalized morphology/correlation QC.
ica_use_ecg_detection = False  # No ECG channel in these recordings.
random_state = 2026
# reject = None  # Decide post-ICA epoch rejection after inspecting first results.
run_source_estimation = False
decode = False
