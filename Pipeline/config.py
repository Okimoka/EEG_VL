# Default settings for data processing and analysis.

from collections.abc import Callable, Sequence
import sys
import os
from typing import Annotated, Any, Literal


from annotated_types import Ge, Interval, Len, MinLen
from mne import Covariance
from mne_bids import BIDSPath

from mne_bids_pipeline.typing import (
     ArbitraryContrast,
     DigMontageType,
     FloatArrayLike,
     PathLike,
 )

import argparse
import numpy as np

bids_root = "./data_split/"

interactive = False


args = [arg for arg in sys.argv if arg.startswith("--task") or not arg.startswith("-")]
parser = argparse.ArgumentParser()
parser.add_argument("ignored", nargs="*")
parser.add_argument(
    "--task", choices=("Gambling","Oddball","Axon"), required=True
)
task = parser.parse_args(args).task 
sessions = "all"

runs = "all" # Always 01 anyway
#crop_runs: tuple[float, float] | None = None
subjects = "all"

ch_types = ["eeg"] # not considering eog channels 
data_type = "eeg" # redundant


eog_channels = None # can be none, it will simply use all actual EOG channels

eeg_reference = "average" # default setting

#TODO check if this actually applies montage
eeg_template_montage = None # custom montage

#drop_channels: Sequence[str] = []

# Read full_config.py for more context
# Might make sense to restrict this?
analyze_channels = "ch_types"


read_raw_bids_verbose = 'error' # Don't log everyting, since files are not perfectly BIDS compliant

#FFT Power spectrum
plot_psd_for_runs = "all" # Default value



######################### Preprocessing #########################


# Read full_config.py for more context
# Paper doesn't mention breaks, but could have still happened since gameplay was quite long
find_breaks: bool = True # Default is false

min_break_duration: float = 15.0 # Default value

# Makes sense, also considering our epoching
t_break_annot_start_after_previous_event: float = 5.0 # Default value
t_break_annot_stop_before_next_event: float = 5.0 # Default value




# ## Bad channel detection
#
# Read full_config.py for context of the used algorithm

# Original paper simply removed most ventral electrodes (leaving 60).
# Is there something like this for eeg channels?
#find_flat_channels_meg: bool = True # Default is False
#find_noisy_channels_meg: bool = True # Default is False


# %%
# ## Maxwell filter

# Read full_config.py for context
# Might make sense to enable? There are many parameters than can be set
# Chose to keep disabled since it doesn't seem to be common? (esp. for EEG data)
# According to the samples: https://mne.tools/mne-bids-pipeline/stable/examples/examples.html
# use_maxwell_filter: bool = False
# mf_st_duration: float | None = None


# ## Filtering & resampling

# ### Filtering
#
# Read full_config.py for context (!!!)
# mne-bids-pipeline does not allow ICA if l_freq is below 1 Hz!! See ica_l_freq in full_config.py
l_freq: float | None = 1 # Paper uses 0.1 Hz, full_config.py recommends 1 Hz
h_freq: float | None = 40 # Paper uses 20 Hz, full_config.py recommends 40 or 120 depending on type of analysis

#TODO Why does this default to None? Is powerline noise filtered elsewhere?
notch_freq: float | Sequence[float] | None = 60 # Default is None
# notch_trans_bandwidth: float = 1.0



# ### Resampling
#
# Not necessary since already 500Hz
# raw_resample_sfreq: float | None = None


# ## Epoching

# Could be useful if event names are not as desired
# rename_events: dict = dict()
# on_rename_missing_events: Literal["ignore", "warn", "raise"] = "raise"

#1) missile launch button press (SHOOT_BUTTON)
#2) collect star (COLLECT_STAR)
#3) collect ammo box (COLLECT_AMMO)
#4) crash into wall (PLAYER_CRASH_WALL)
#5) crash into enemy (PLAYER_CRASH_ENEMY)
#6) missile hit enemy (MISSILE_HIT_ENEMY)

"""
The paper describes issues with the "PLAYER_CRASH_WALL" and "COLLECT_STAR" events often occuring in rapid succession, where only the first one is kept.
The filtering has already been done in step 01, and the offending events are now marked as "IGNORE PLAYER_CRASH_WALL" and "IGNORE COLLECT_STAR"
"""

conditions = []

if(task == "Gambling"): #"Wins and losses" Exemplar events
    conditions = ["GAMBLING LOSS", "GAMBLING WIN"]
elif(task == "Oddball"): #"Target detection" Exemplar events
    conditions = ["ODDBALL STANDARD", "ODDBALL RARE"]
elif(task == "Axon"):
    #There are also "GAME OVER" and "GAME START", but they occur only once (?)
    conditions = ["SHOOT_BUTTON", "MISSILE_HIT_ENEMY", "COLLECT_STAR", "COLLECT_AMMO", "PLAYER_CRASH_ENEMY", "PLAYER_CRASH_WALL"]

#TODO Check whether these epoch timings from the paper make sense. They seem unusually long?
epochs_tmin: float = -2
epochs_tmax: float = 2

#Also from the paper
baseline: tuple[float | None, float | None] | None = (-0.2, 0)


# ## Artifact removal

# ### SSP, ICA, and artifact regression

#TODO check if this performs the intended behavior (removing eye artifacts from eeg by using the eog channels)
# For some reason, this throws an error :(
#regress_artifact = {
#    "picks": "eeg", 
#    "picks_artifact": ["HEOG", "VEOG"]
#}

#TODO maybe SSP could make sense too?
# For SSP there are some more variables aswell
spatial_filter: Literal["ssp", "ica"] | None = "ica"
# """
# !!! warning "ICA requires manual intervention!"
#     After the automatic ICA component detection step, review each subject's
#     `*_report.html` report file check if the set of ICA components to be removed
#     is correct. Adjustments should be made to the `*_proc-ica_components.tsv`
#     file, which will then be used in the step that is applied during ICA.

#     ICA component order can be considered arbitrary, so any time the ICA is
#     re-fit – i.e., if you change any parameters that affect steps prior to
#     ICA fitting – this file will need to be updated!
# """

#TODO check percentage of rejected epochs. Paper states 3.27%
ica_reject: dict[str, float] | Literal["autoreject_local"] | None = "autoreject_local"


# Can go as high as 3000 for convergence
# ica_max_iterations: int = 500

# Read full_config.py for context
# None == 0.9999 might be optimal if only artifact removal is of interest
# ica_n_components: float | int | None = None


#TODO Browse raw data to find actual eog and ecg thresholds
#Read full_config.py for tips on how to do this
# ica_ecg_threshold: float = 0.1
# """
# The cross-trial phase statistics (CTPS) threshold parameter used for detecting
# ECG-related ICs.
# """

# ica_eog_threshold: float = 3.0
# """
# The threshold to use during automated EOG classification. Lower values mean
# that more ICs will be identified as EOG-related. If too low, the
# false-alarm rate increases dramatically.
# """

#For rejection after ICA
#Probably makes the most sense for data with moderate amounts of artifacts
reject: dict[str, float] | Literal["autoreject_global", "autoreject_local"] | None = "autoreject_local"#(
#     None
# )


# In the paper a median of 2 electrodes (range 1 to 3) was interpolated, so add these values
autoreject_n_interpolate: FloatArrayLike = [1,2,3, 4, 8, 16] # Default [4, 8, 16]
# """
# The maximum number of bad channels in an epoch that `autoreject` local will try to
# interpolate. The optimal number among this list will be estimated using a
# cross-validation procedure; this means that the more elements are provided here, the
# longer the `autoreject` run will take. If the number of bad channels in an epoch
# exceeds this value, the channels won't be interpolated and the epoch will be dropped.

# Read full_config.py for more context
# Has many more parameters
decode = False # Default is True

# ## Time-frequency analysis

# Just put in all events for now (?)
time_frequency_conditions: Sequence[str] = conditions

#TODO look at the data to find optimal values
# time_frequency_freq_min: float | None = 8
# time_frequency_freq_max: float | None = 40

time_frequency_freq_min = 1  # minimum frequency in Hz
time_frequency_freq_max = 50  # maximum frequency in Hz
#1 to 50 Hz in 49 logarithmically spaced steps
frequencies = np.logspace(np.log10(time_frequency_freq_min), np.log10(time_frequency_freq_max), 49)

# sigma
time_frequency_cycles = 4 / (2 * np.pi * frequencies)

# time_frequency_subtract_evoked: bool = False

time_frequency_baseline: tuple[float, float] | None = (-0.3, -0.2)
time_frequency_baseline_mode: str = "logratio" # Default is "mean"
time_frequency_crop: dict | None = dict(tmin=-0.5, tmax=1, fmin=1, fmax=50) # Default in None



n_jobs: int = os.cpu_count()
# """
# Specifies how many subjects you want to process in parallel. If `1`, disables
# parallel processing.
# """
