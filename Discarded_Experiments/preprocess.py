
"""


Attempts:
- Remove bad channels by analyzing kurtosis of channels
- Splitting and proper reading of data, now done in 01_split_data.py



Data Sampled across .01 to 100Hz
Sampling rate 500Hz
Reference CPz
64 Channels
40 ms delay from trigger onset to visual change on monitor

VEOG und HEOG (vertical and horizontal electrooculogram) are the two EOG channels meaning they are used to measure eye movements
and they are used to remove eye movement artifacts from the EEG data
"""

import matplotlib.pyplot as plt
import time
import sys
import os
import logging
import ccs_eeg_utils
import pandas as pd
import numpy as np
import mne, mne_bids, ccs_eeg_utils
from scipy.stats import kurtosis

#Load the data
from mne_bids import (BIDSPath, read_raw_bids)


bids_root = "./data/"

#none of the channels are marked as bad in the dataset (raw.info["bads"])



def get_bad_channels(raw):
    # Get EEG channel data
    eeg_data = raw.get_data(picks='eeg')

    # Compute kurtosis for each channel
    kurt_values = kurtosis(eeg_data, axis=1)

    # Create a list of channels with their corresponding kurtosis values
    channels_kurtosis = list(zip(raw.ch_names, kurt_values))

    # Sort the list by kurtosis values in descending order
    sorted_channels_kurtosis = sorted(channels_kurtosis, key=lambda x: abs(x[1]), reverse=True)

    # Convert to DataFrame for better display
    df_kurtosis = pd.DataFrame(sorted_channels_kurtosis, columns=['Channel', 'Kurtosis'])

    # Display the sorted list
    print(df_kurtosis)

    # Plot distributions of all channels
    plt.figure(figsize=(15, 8))

    # Plot each channel's data
    for idx, channel_data in enumerate(eeg_data):
        plt.plot(channel_data, label=f'{raw.ch_names[idx]}')

    plt.title("EEG Channel Distributions")
    plt.xlabel("Time Points")
    plt.ylabel("Amplitude (µV)")
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1), ncol=2)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()



#bids_path = BIDSPath(subject="001",task="Gambling",run="01", datatype='eeg', suffix='eeg', root=bids_root)
#bids_path = BIDSPath(subject="001",task="Oddball",run="01", datatype='eeg', suffix='eeg', root=bids_root)
#bids_path = BIDSPath(subject="001",task="Axon",run="01", datatype='eeg', suffix='eeg', root=bids_root)


for x in range(1,17): #actually (1,17)
    with mne.utils.use_log_level('error'):
        subject_id = '0'+str(x) if x > 9 else '00'+str(x)
        print("Subject: ", subject_id)

        # data/sub-001/eeg/sub-001_task-ContinuousVideoGamePlay_run-01_eeg.set
        bids_path = BIDSPath(subject=subject_id,task="ContinuousVideoGamePlay",run="01",
                            datatype='eeg', suffix='eeg',
                            root=bids_root)
        #read the file
        raw = read_raw_bids(bids_path)
        raw.load_data()

        # channels are all n/a in the dataset, but all should be eeg (except for VEOG and HEOG channels which are eog)
        channel_types = {ch: 'eeg' if not ch.endswith("EOG") else 'eog' for ch in raw.ch_names}
        raw.set_channel_types(channel_types)


        # Coordinates in _electrodes.tsv are not recognized automatically (?), so they are read manually like this
        # Positions of the eog channels are automatically ignored
        electrodes_path = bids_path.copy().update(suffix='electrodes', extension='.tsv')
        electrodes_df = pd.read_csv(electrodes_path.fpath, sep='\t')
        ch_pos = {row['name']: (row['x'], row['y'], row['z']) for _, row in electrodes_df.iterrows()}
        custom_montage = mne.channels.make_dig_montage(ch_pos, coord_frame='head')
        raw.set_montage(custom_montage, match_case=False)
        #mne.viz.plot_montage(custom_montage)
        #raw.plot_sensors(ch_type="eeg", block=True)



    annotations = raw.annotations 

    # Find all events with "STATUS" in the description
    # Ignore them if they are the last or first event
    status_events = [
        annotation for annotation in annotations
        if "STATUS" in annotation['description'] and annotation != annotations[0] and annotation != annotations[-1]
    ]

    status_times = [annotation['onset'] for annotation in status_events]

    if(len(status_events) == 0):
        # Oddball task and gambling task are split by "STATUS" event
        # If there is only one "STATUS" event, the subject did not complete the gambling task
        print("Skipping subject who did not complete the gambling task")
        continue

    #TODO For now just skip difficult subjects
    if(subject_id == "008"):
        # For subject 008, part of the Axon gameplay leaked into Oddgamble events (separated by "STATUS")
        continue

    if(subject_id == "006"):
        # Subject 006 played 3 sessions of gambling. Maybe just merge them?
        continue


    # Split the raw data
    raw_oddball = raw.copy().crop(tmin=raw.times[0], tmax=status_times[0])
    raw_gambling = raw.copy().crop(tmin=status_times[0], tmax=raw.times[-1])

    #"EEG data were re-referenced to an average reference and CPz was re-created via the EEGLab pop_reref function"
    # TODO does this make sense? CPz does not have a coordinate in the montage
    raw_oddball, _ = mne.set_eeg_reference(raw_oddball, ref_channels='average')
    raw_oddball = mne.add_reference_channels(raw_oddball, ref_channels=['CPz'])


    #bad_channels_oddball = get_bad_channels(raw_oddball)





    ## Get channel locations
    #montage = raw_oddball.get_montage()

    ## Extract positions of all channels
    #ch_positions = montage.get_positions()['ch_pos']

    ## Separate good and bad channels
    #good_positions = {ch: pos for ch, pos in ch_positions.items() if ch not in raw_oddball.info['bads']}
    #bad_positions = {ch: pos for ch, pos in ch_positions.items() if ch in raw_oddball.info['bads']}

    ## Visualize using MNE's built-in plotting function
    #fig = plt.figure(figsize=(10, 8))
    #mne.viz.plot_montage(montage, show_names=True, show=True)

    ## Highlight bad channels manually if desired
    #if bad_positions:
    #    bad_montage = mne.channels.make_dig_montage(bad_positions, coord_frame='head')
    #    mne.viz.plot_montage(bad_montage, show_names=True, sphere=0.1, color='red')




    #print(len(blocks))


    # Print the result
    #print("Indexes of events with 'STATUS':", status_event_indexes)
    #print("--------------------------")




    #print(events2)
    #boundary_idx = [e[0] for e in events if event_id['STATUS/boundary'] == e[2]]
    #print("Boundary index: ", boundary_idx)
    
    #events, event_id = mne.events_from_annotations(raw)
    #print(event_id)  # See what event types exist
    #print(events)    # Inspect all detected events
    
    #print(list(raw.annotations))
    """
    blocks = prep.split_in_blocks(raw.copy())

    if os.path.isfile("./data/bad_channels.json"):
        bads = json.load(open("./data/bad_channels.json"))
        blocks = prep.set_bad_channels_from_json(blocks, bads)
    else:
        bads = misc.create_bad_json_structure()

    if os.path.isfile("./data/bad_ica_components.json"):
        ica_bads = json.load(open("./data/bad_ica_components.json"))
    else:
        ica_bads = misc.create_bad_json_structure()
    
    # deprecated, do not use
    if Z_SCORE_REJECT == True:
        for b in blocks:
            # reject by z-score (autoreject is more sophisticated, but only works on epochs and is really slow)
            b.info['bads'].extend(prep.mark_bad_channels_by_z_score(b, threshold=5.0, window_size=10000))
            print(f"Bad channels: {b.info['bads']}")

    if PYPREP_REJECT == True:
        for b in blocks:
            nc = NoisyChannels(b, random_state=42) # of course it has to be 42!
            nc.find_all_bads(ransac=False)
            b.info['bads'].extend(nc.get_bads())
            print(f"Bad channels: {b.info['bads']}")

    if PROMPT_BADS == True:
        for b in blocks:
            b.plot(n_channels=64)
            plt.show(block=True)
            bads[f"sub-{SUBJECT}"][f"{blocks.index(b)+1}"] = b.info['bads']
    
    
    with open("./data/bad_channels.json", "w") as f:
        json.dump(bads, f)

    ica_blocks = []

    # separate preprocessing in 8 blocks as boundary events occured 8 times in the whole recording
    for b in blocks:
        b.interpolate_bads()
        prep_block = prep.basic_preprocessing(b.copy())

        # ICA
        ica_block = prep.get_ica(prep_block, ica_bads, f"{blocks.index(b)+1}", montage)
        ica_blocks.append(ica_block)
    
    with open("./data/bad_ica_components.json", "w") as f:
        json.dump(ica_bads, f)

    prep_raw = mne.concatenate_raws(ica_blocks)
    misc.save_preprocessed_data(f"./data/fifs/processed_{SUBJECT}_raw.fif", prep_raw)

    else:
    prep_raw = mne.io.read_raw_fif(f"./data/fifs/processed_{SUBJECT}_raw.fif", preload=True)

    prep_raw.info
    """




#ccs_eeg_utils.read_annotations_core(bids_path,raw)
#annotations = mne.read_annotations("sub-{}_task-P3_badannotations.csv".format(subject_id))
#raw.annotations.append(annotations.onset,annotations.duration,annotations.description)


#raw.load_data()
#raw.filter(0.5,50, fir_design='firwin')
#raw.plot(n_channels=len(raw.ch_names), block=True)#,scalings =40e-6)


# See below










#print(raw.info['chs'])



#ccs_eeg_utils.read_annotations_core(bids_path,raw)
#raw.load_data()