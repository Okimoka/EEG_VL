import mne
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np

def analyze_subjects(subject_ids, task_choice, electrode):
    grand_avg_data_standard = []
    grand_avg_data_oddball = []
    
    for subject_id in subject_ids:
        if task_choice == 'Gambling':
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-Gambling_proc-ica_epo.fif'
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-Gambling_proc-clean_epo.fif'
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/bids_output/sub-{subject_id}/eeg/sub-{subject_id}_task-ContinuousVideoGamePlay_run-1_eeg.set'
            filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-OddGamble_proc-clean_epo.fif'
            standard_label = 'GAMBLING LOSS'
            oddball_label = 'GAMBLING WIN'
        elif task_choice == 'Oddball':
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-Oddball_proc-ica_epo.fif'
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-Oddball_proc-clean_epo.fif'
            #filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/bids_output/sub-{subject_id}/eeg/sub-{subject_id}_task-ContinuousVideoGamePlay_run-1_eeg.set'
            filepath = f'/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/sub-{subject_id}/eeg/sub-{subject_id}_task-OddGamble_proc-clean_epo.fif'
            standard_label = 'ODDBALL STANDARD'
            oddball_label = 'ODDBALL RARE'
        else:
            print("Invalid task choice. Please choose 'Gambling' or 'Oddball'.")
            return

        epochs = None
        if(filepath.endswith(".set")):

            # Load the raw data
            raw = mne.io.read_raw_eeglab(filepath, preload=True)
            
            # Find events
            standard_event_id = 2
            oddball_event_id = 3
            events, event_id = mne.events_from_annotations(raw)
            if task_choice == 'Gambling':
                standard_event_id = 6
                oddball_event_id = 7

            print(events)
            print(event_id)
            # Specify the event IDs that correspond to standard and oddball labels
            #standard_event_id = event_id[standard_label]
            #oddball_event_id = event_id[oddball_label]

            # Create epochs
            epochs = mne.Epochs(raw, events, event_id={standard_label: standard_event_id, oddball_label: oddball_event_id}, tmin=-2, tmax=2, preload=True)
        else:
            epochs = mne.read_epochs(filepath, preload=True)

        # Find index of the electrode
        ch_idx = epochs.ch_names.index(electrode)

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
        smooth_standard = grand_avg_standard
        smooth_oddball = grand_avg_oddball


        times = epochs.times  # Time vector can be assumed to be the same for all

        # Plotting for each subject
        #plt.figure(figsize=(10, 6))
        #plt.plot(times, smooth_standard * 1e6, label=f'{standard_label} at {electrode} (Smoothed)', color='blue')
        #plt.plot(times, smooth_oddball * 1e6, label=f'{oddball_label} at {electrode} (Smoothed)', color='red')
        #plt.xlabel('Time (s)')
        #plt.ylabel('Amplitude (μV)')
        #plt.title(f'{standard_label} vs {oddball_label} at {electrode} for Subject {subject_id}')
        #plt.legend()
        #plt.grid(True)
        #plt.show()

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
    ax = plt.gca()
    ax.set_ylim([-5, 10])
    plt.show()



# Example usage:
#subject_ids = ["001","004","007","009","010",'011','012','013','014','015','016','017']#input("Enter the subject IDs separated by space: ").split()
subject_ids = [str(i).zfill(3) for i in range(1, 18)]
task_choice = "Oddball"
electrode = "Pz"
analyze_subjects(subject_ids, task_choice, electrode)
