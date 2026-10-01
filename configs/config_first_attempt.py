from pathlib import Path

_here = Path(__file__).resolve().parents[1]
bids_root = _here / "prepared_bids"
deriv_root = _here / "logs" / "first_attempt"
subjects = "all"
runs = ["01"]
task = "ContinuousVideoGamePlay"
ch_types = ["eeg"]
conditions = ['ODDBALL STANDARD', 'ODDBALL RARE']
l_freq = 0.1
h_freq = 40
random_state = 2026
epochs_tmin = -0.5
epochs_tmax = 1.0
run_source_estimation = False
decode = False
