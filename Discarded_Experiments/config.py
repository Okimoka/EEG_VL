"""
I think this is a (modified) example configuration file from the MNE-BIDS-pipeline examples
https://mne.tools/mne-bids-pipeline/stable/examples/examples.html
"""


import argparse
import sys

import mne

bids_root = "./data_split"
deriv_root = "./results"

#/sub-001/eeg/sub-001_task-P3_run-01

# Find the --task option
args = [arg for arg in sys.argv if arg.startswith("--task") or not arg.startswith("-")]
parser = argparse.ArgumentParser()
parser.add_argument("ignored", nargs="*")
parser.add_argument(
    "--task", choices=("N400", "ERN", "LRP", "MMN", "N2pc", "N170", "Oddball"), required=True
)
task = parser.parse_args(args).task

conditions = ["STATUS"]
#contrasts = [("STATUS")]


subjects = ["001"]

ch_types = ["eeg"]
interactive = False

raw_resample_sfreq = 128
# Suppress "Data file name in EEG.data (sub-019_task-ERN_eeg.fdt) is incorrect..."
read_raw_bids_verbose = "error"

eeg_template_montage = mne.channels.make_standard_montage("standard_1005")
#eeg_bipolar_channels = {
#    "HEOG": ("HEOG_left", "HEOG_right"),
#    "VEOG": ("VEOG_lower", "FP2"),
#}
#drop_channels = ["HEOG_left", "HEOG_right", "VEOG_lower"]
#eog_channels = ["HEOG", "VEOG"]

l_freq = 0.1
h_freq = None
notch_freq = 60


if task == "N400":  # test autoreject local without ICA
    spatial_filter = None
    reject = "autoreject_local"
    autoreject_n_interpolate = [2, 4]
elif task == "N170":  # test autoreject local before ICA
    spatial_filter = "ica"
    ica_reject = "autoreject_local"
    reject = "autoreject_global"
    autoreject_n_interpolate = [2, 4]
else:
    spatial_filter = "ica"
    ica_reject = dict(eeg=350e-6, eog=500e-6)
    reject = "autoreject_global"

# These settings are only used for the cases where spatial_filter="ica"
ica_max_iterations = 1000
ica_eog_threshold = 2
ica_decim = 2  # speed up ICA fitting

run_source_estimation = False
on_rename_missing_events = "ignore"

n_jobs = 4


eeg_reference = ["CPz"]