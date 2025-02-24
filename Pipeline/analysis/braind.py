import os
import pickle
import sys
import mne
from braindecode import EEGClassifier
from skorch.callbacks import EarlyStopping
from torch.optim import Adam
import torch
from sklearn.model_selection import train_test_split
import numpy as np

from braindecode.datasets import (
    create_from_mne_raw, create_from_mne_epochs)

from braindecode.models import ShallowFBCSPNet, EEGNetv4, EEGInceptionERP, Deep4Net, TIDNet
from sklearn.model_selection import StratifiedKFold

task = "Oddball"
home_path = "/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/data_split/derivatives/"
cache_path = "/home/okk/Documents/Letztes_Semester/EEG_VL/EVERYTHING_IMPORTANT/Pipeline/analysis/classifier_cache/"

os.makedirs(cache_path, exist_ok=True)

subjects = [f"{i:03}" for i in range(1, 18)]

for subject in subjects:

    #if(task == "Gambling"):
    #    #these subjects did not gamble
    #    if(subject == "002" or subject == "003" or subject == "005"):
    #        continue

    # Count occurrences of each event


    print("SUBJECT: ", subject)
    game_file_path = home_path + f"sub-{subject}/eeg/sub-{subject}_task-Axon_proc-clean_epo.fif"

    #only process subjects that have the game file
    if not os.path.exists(game_file_path):
        print(f"Skipping subject {subject} due to missing game file.")
        continue

    #game_file_path = 'sub-001_task-Axon_proc-ica_epo.fif'
    # Load the data
    file_path = home_path + f"sub-{subject}/eeg/sub-{subject}_task-OddGamble_proc-clean_epo.fif"

    cache_filename = f"{task.lower()}_classifier_sub-{subject}.pkl"
    full_cache_path = os.path.join(cache_path, cache_filename)
    #file_path = home_path + f"sub-{subject}/eeg/sub-{subject}_task-Gambling_proc-clean_epo.fif"
    #file_path = f'sub-{subject}_task-Gambling_proc-clean_epo.fif'
    #if(task == "Oddball"):
    #    file_path = f'sub-{subject}_task-Oddball_proc-clean_epo.fif'
    #axon_path = 'sub-001_task-Axon_proc-ica_epo.fif'

    #gambling_epochs = mne.read_epochs(gambling_path)
    clf = None


    epochs = mne.read_epochs(file_path)
    game_epochs = mne.read_epochs(game_file_path)

    print("Game Event IDs:", game_epochs.event_id)

    event_counts = {event_id: np.sum(epochs.events[:, -1] == event_id) for event_id in [1,2,4,5]}
    # Check if any event has fewer than 5 occurrences
    if any(count < 10 for count in event_counts.values()):
        print(f"Skipping subject {subject} due to insufficient events.")
        continue  # Skip to the next subject


    # Create a mapping for original event IDs to binary labels
    event_dict = {'GAMBLING WIN': 2, 'GAMBLING LOSS': 1}
    if(task == "Oddball"):
        event_dict = {'ODDBALL RARE': 4, 'ODDBALL STANDARD': 5}

    epochs = epochs[np.isin(epochs.events[:, -1], list(event_dict.values()))]


    # Define correct mapping to convert event IDs
    event_id_mapping = {1: 0, 2: 1}
    if(task == "Oddball"):
        event_id_mapping = {4: 0, 5: 1}

    # Convert the labels
    binary_labels = np.array([event_id_mapping[e[-1]] for e in epochs.events])
    print(binary_labels)
    #sys.exit(0)
    #print(np.array([event_id_mapping[e[-1]] for e in epochs.events]))
    #print(np.array([event_id_mapping[e] for e in epochs.events[:, -1]]))
    # Overwrite the event IDs in the epochs object
    epochs.events[:, -1] = np.array([event_id_mapping[e] for e in epochs.events[:, -1]])

    # Overwrite the event_id dictionary to reflect the new mapping
    epochs.event_id = {'GAMBLING WIN': 1, 'GAMBLING LOSS': 0}
    if(task == "Oddball"):
        epochs.event_id = {'ODDBALL RARE': 0, 'ODDBALL STANDARD': 1}

    #list_of_epochs = [mne.Epochs(raw, [[0, 0, 0]], tmin=0, baseline=None) for raw in parts]

    windowsize = 3
    samplingrate = 500

    window_ds = create_from_mne_epochs(
        [epochs],
        window_size_samples=windowsize*samplingrate + 1,  # adjust based on your sampling rate and desired window size
        window_stride_samples=windowsize*samplingrate + 1,  # adjust stride to control overlap
        drop_last_window=False
    )
    #print(epochs.get_data().shape)
    #print(f"Total windows: {len(window_ds)}")
    #print("Example window shape:", window_ds[0][0].shape)
    #print(f"Number of samples in dataset: {len(window_ds)}")
    #print(f"Number of labels: {len(binary_labels)}")

    print("Event IDs:", epochs.event_id)
    print("Event values:", np.unique(epochs.events[:, -1]))
    #exit(0)

    # Define the number of EEG channels and the number of time samples in your input window
    n_channels = window_ds[0][0].shape[0]  # Number of EEG channels
    n_times = window_ds[0][0].shape[1]     # Number of time samples

    # Initialize the ShallowFBCSPNet model
    model = ShallowFBCSPNet(
        n_chans=n_channels,
        n_outputs=2,  # Assuming binary classification
        n_times=n_times,
        final_conv_length='auto'
    )
    clf = EEGClassifier(
        model,
        criterion=torch.nn.CrossEntropyLoss,
        optimizer=Adam,
        train_split=None,
        optimizer__lr=0.01,
        batch_size=64,
        device='cpu',
        verbose=False
    )


    result = []
    # Implement StratifiedKFold
    skf = StratifiedKFold(n_splits=5)
    for fold, (train_idx, test_idx) in enumerate(skf.split(window_ds, binary_labels)):

        cache_filename = f"{task.lower()}_classifier_sub-{subject}_fold-{fold+1}.pkl"
        full_cache_path = os.path.join(cache_path, cache_filename)

        train_dataset = torch.utils.data.Subset(window_ds, train_idx)
        test_dataset = torch.utils.data.Subset(window_ds, test_idx)

        if os.path.exists(full_cache_path):
            print(f"Loading cached classifier for {subject}, fold {fold+1}")
            with open(full_cache_path, 'rb') as f:
                clf = pickle.load(f)
        else:
            # Initialize and train the classifier if no cache is available

            # Assume clf has been initialized before the fold loop if it's not specific to each fold
            clf.fit(train_dataset, y=binary_labels[train_idx])

            # Cache the trained classifier for this fold
            with open(full_cache_path, 'wb') as f:
                pickle.dump(clf, f)
            print(f"Cached classifier for {subject}, fold {fold+1}")

        # Evaluate the model using test data
        accuracy = clf.score(test_dataset, y=binary_labels[test_idx])
        result.append(f"Fold {fold+1} - Accuracy: {accuracy * 100:.2f}%")

    #Document results in text file
    with open(f"results/results_{subject}_{task}.txt", "a") as f:
        f.write(f"Subject {subject}\n")
        f.write(f"Task: {task}\n")
        f.write(f"Event counts: {event_counts}\n")
        f.write(f"Results:\n")
        f.write("\n".join(result))
        f.write("\n\n")
        print("Modified results file")



#game_event_dict = {'MISSILE_HIT': 1, 'CRASH': 0}
#game_labels = np.array([game_event_dict[e[-1]] for e in game_epochs.events])