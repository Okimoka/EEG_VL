"""Broadband continuous ICA application only; no new fit or epoch rejection."""
from pathlib import Path
_root = Path(__file__).resolve().parent
bids_root = _root / 'prepared_native_bids'
deriv_root = _root / 'artifacts' / 'ica_broadband'
subjects = [f'{number:03d}' for number in range(1, 18)]
runs = ['01', '02']
task = 'ContinuousVideoGamePlay'
ch_types = ['eeg']
data_type = 'eeg'
n_jobs = 4
mne_log_level = 'info'
l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0
# raw_resample_sfreq = None  # Preserve 500 Hz, including EOG.
# eeg_reference = 'average'  # Same good-channel reference as the saved ICA.
# Required event configuration for config validation; no epochs are created here.
conditions = ['ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING WIN', 'GAMBLING LOSS',
              'PLAYER_CRASH_WALL', 'COLLECT_STAR', 'MISSILE_HIT_ENEMY',
              'PLAYER_CRASH_ENEMY', 'COLLECT_AMMO', 'SHOOT_BUTTON']
baseline = None
spatial_filter = 'ica'
ica_use_icalabel = True  # Native apply uses the fitted average-reference convention.
ica_use_eog_detection = False
ica_use_ecg_detection = False
ica_algorithm = 'extended_infomax'
ica_h_freq = 100.0  # Describes the saved model; no ICA fit in this branch.
ica_reject = 'autoreject_local'  # Describes the original training selection only.
reject = None  # No analysis epochs or AutoReject in the continuous carrier.
random_state = 2026
run_source_estimation = False
decode = False
