from pathlib import Path
_root = Path(__file__).resolve().parents[1]
bids_root = _root / 'reviewed_bids'
deriv_root = _root / 'logs' / 'ica_pilot'
subjects = ['001', '007']
runs = ['01', '02']
task = 'ContinuousVideoGamePlay'
ch_types = ['eeg']
data_type = 'eeg'
n_jobs = 4
memory_location = False
l_freq = 0.1
h_freq = 100.0
notch_freq = 60.0
conditions = ['ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING WIN', 'GAMBLING LOSS',
              'PLAYER_CRASH_WALL', 'COLLECT_STAR', 'MISSILE_HIT_ENEMY',
              'PLAYER_CRASH_ENEMY', 'COLLECT_AMMO', 'SHOOT_BUTTON']
epochs_tmin = -2.0
epochs_tmax = 2.0
baseline = None
spatial_filter = 'ica'
ica_reject = 'autoreject_local'
ica_algorithm = 'extended_infomax'
ica_h_freq = 100.0
ica_use_icalabel = True
ica_use_eog_detection = False
ica_use_ecg_detection = False
random_state = 2026
run_source_estimation = False
decode = False
