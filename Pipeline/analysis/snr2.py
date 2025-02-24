import mne
from scipy.io import loadmat
import numpy as np

def load_eeg_data(filepath):
    # Load EEG data using MNE (assuming EEGLAB .set format)
    raw = mne.io.read_raw_eeglab(filepath, preload=True)
    return raw

def calculate_pnr(raw, cleaned, event_id, tmin, tmax, baseline):
    # Find events in the raw data
    events, event_ids = mne.events_from_annotations(raw)
    # Define epochs around the event of interest
    epochs = mne.Epochs(raw, events, event_id=event_id, tmin=tmin, tmax=tmax,
                        baseline=baseline, preload=True)
    cleaned_epochs = mne.Epochs(cleaned, events, event_id=event_id, tmin=tmin, tmax=tmax,
                                baseline=baseline, preload=True)
    
    # Calculate ERP from the cleaned data
    evoked_cleaned = cleaned_epochs.average()
    # Peak amplitude in the signal
    peak_amp = np.max(np.abs(evoked_cleaned.data))
    
    # Noise level from the baseline period in the raw data
    baseline_epochs = mne.Epochs(raw, events, event_id=event_id, tmin=baseline[0], tmax=baseline[1],
                                 baseline=baseline, preload=True)
    noise_level = np.std(baseline_epochs.get_data())

    # Calculate Peak-to-Noise Ratio
    pnr = peak_amp / noise_level
    return pnr

# File paths
raw_filepath = 'raw.set'
cleaned_filepath = 'raw2.set'

# Load data
raw_data = load_eeg_data(raw_filepath)
cleaned_data = load_eeg_data(cleaned_filepath)

# Parameters for ODDBALL RARE
event_id = {'ODDBALL RARE': 3}  # using the ID for rare oddball events
tmin, tmax = -0.2, 1  # time window around the event (200 ms pre-stimulus to 800 ms post-stimulus)
baseline = (-0.2, 0)  # baseline period (200 ms pre-stimulus)

# Calculate SNR
snr = calculate_pnr(raw_data, cleaned_data, event_id, tmin, tmax, baseline)
print(f"Peak-to-Noise Ratio (PNR) is: {snr:.2f}")