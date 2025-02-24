import mne
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np
import pandas as pd

from mne.preprocessing import ICA
from mne_icalabel import label_components


def analyze_subjects(subject_ids, task_choice, electrode):
    grand_avg_data_standard = []
    grand_avg_data_oddball = []

    for subject_id in subject_ids:
        filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data/sub-{subject_id}/eeg/sub-{subject_id}_task-ContinuousVideoGamePlay_run-01_eeg.set'  # Update the path accordingly
        electrodes_path = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data/sub-{subject_id}/eeg/sub-{subject_id}_task-ContinuousVideoGamePlay_run-01_electrodes.tsv'  # Update the path accordingly
        # Load raw data                                                                                                         sub-001_task-ContinuousVideoGamePlay_run-01_electrodes.tsv
        raw = mne.io.read_raw_eeglab(filepath, preload=True)


        #electrodes_path = bids_path.copy().update(suffix='electrodes', extension='.tsv')
        electrodes_df = pd.read_csv(electrodes_path, sep='\t')
        ch_pos = {row['name']: (row['x'], row['y'], row['z']) for _, row in electrodes_df.iterrows()}
        custom_montage = mne.channels.make_dig_montage(ch_pos, coord_frame='head')
        raw.set_montage(custom_montage, match_case=False)

        # Preprocessing: Filtering and re-referencing
        raw.filter(l_freq=1.0, h_freq=100.0)  # Band-pass filter
        raw.set_eeg_reference(ref_channels="average")  # Average reference
        # Run ICA for artifact removal
        ica = ICA(n_components=15, random_state=97)
        ica.fit(raw)
        ica_labels = label_components(raw, ica, method="iclabel")
        artifact_inds = [idx for idx, label in enumerate(ica_labels["labels"]) if label not in ["brain", "other"]]
        ica.exclude = artifact_inds
        ica.apply(raw)






        # Find event information
        events, event_id = mne.events_from_annotations(raw)

        # Define epochs around events
        epochs = mne.Epochs(raw, events, event_id=event_id, tmin=-1, tmax=1, preload=True)

        # Find index of the electrode
        ch_idx = epochs.ch_names.index(electrode)

        standard_label = 'S 22'  # Update with your actual event label
        oddball_label = 'S 25'   # Update with your actual event label

        standard_epochs = epochs[standard_label].get_data()[:, ch_idx, :]
        oddball_epochs = epochs[oddball_label].get_data()[:, ch_idx, :]

        # Compute the grand average for each subject
        grand_avg_standard = np.mean(standard_epochs, axis=0)
        grand_avg_oddball = np.mean(oddball_epochs, axis=0)
        # Accumulate data for grand average across subjects
        grand_avg_data_standard.append(grand_avg_standard)
        grand_avg_data_oddball.append(grand_avg_oddball)
        # Apply Savitzky-Golay filter for smoothing individual subject data
        #smooth_standard = savgol_filter(grand_avg_standard, window_length=51, polyorder=3)
        #smooth_oddball = savgol_filter(grand_avg_oddball, window_length=51, polyorder=3)
        times = epochs.times  # Time vector can be assumed to be the same for all
        # Plotting for each subject
        ##plt.figure(figsize=(10, 6))
        ##plt.plot(times, grand_avg_standard * 1e6, label=f'{standard_label} at {electrode} (Smoothed)', color='blue')
        ##plt.plot(times, grand_avg_oddball * 1e6, label=f'{oddball_label} at {electrode} (Smoothed)', color='red')
        ##plt.xlabel('Time (s)')
        ##plt.ylabel('Amplitude (μV)')
        ##plt.title(f'{standard_label} vs {oddball_label} at {electrode} for Subject {subject_id}')
        ##plt.legend()
        ##plt.grid(True)
        ##plt.show()

    # Calculate the grand average across all subjects
    grand_avg_standard_all = np.mean(grand_avg_data_standard, axis=0)
    grand_avg_oddball_all = np.mean(grand_avg_data_oddball, axis=0)

    # Apply Savitzky-Golay filter for smoothing grand average data
    #smooth_standard_all = savgol_filter(grand_avg_standard_all, window_length=51, polyorder=3)
    #smooth_oddball_all = savgol_filter(grand_avg_oddball_all, window_length=51, polyorder=3)

    # Plotting grand average across all subjects
    plt.figure(figsize=(10, 6))
    plt.plot(times, grand_avg_standard_all * 1e6, label=f'Grand Average {standard_label} at {electrode} (Smoothed)', color='blue')
    plt.plot(times, grand_avg_oddball_all * 1e6, label=f'Grand Average {oddball_label} at {electrode} (Smoothed)', color='red')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude (μV)')
    plt.title(f'Grand Average {standard_label} vs {oddball_label} at {electrode} Across All Subjects')
    plt.legend()
    plt.grid(True)
    plt.show()
    # Compute and plot averages as in your existing code
    # Rest of the code remains similar to your original script, with necessary updates

# Example usage
subject_ids = ["001","004","007","009","010",'011','012','013','014','015','016','017']
task_choice = "Gambling"
electrode = "Cz"
analyze_subjects(subject_ids, task_choice, electrode)