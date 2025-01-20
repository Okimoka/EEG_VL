
"""
Currently, Oddball and Gamble datasets are split by the "STATUS" event in the annotations.
This script separates them into two tasks
"""

import matplotlib.pyplot as plt
import time
import sys
import os
import json
import logging
import ccs_eeg_utils
import pandas as pd
import numpy as np
import mne, mne_bids, ccs_eeg_utils
from scipy.stats import kurtosis


#Load the data
from mne_bids import (BIDSPath, read_raw_bids, write_raw_bids)

# Helper method to update the .json file of the subject tasks.
# When processing the data in an mne compliant format for splitting etc., the json files
# No longer retain their original information
def update_json(bids_path, eeg_placement_scheme, eeg_ground, eeg_reference):
    json_path = bids_path.copy().update(extension='.json').fpath
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            data = json.load(f)
        # Update the fields
        data['EEGPlacementScheme'] = eeg_placement_scheme
        data['EEGGround'] = eeg_ground
        data['EEGReference'] = eeg_reference
        
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=4)



bids_root = "./data/"
bids_output_root = "./data_split/"

for x in range(1,18): #for all 17 subjects
    for run in range(1,3): #for both runs
        with mne.utils.use_log_level('error'):
            subject_id = f'{x:03}'
            print("Subject: ", subject_id)

            # data/sub-001/eeg/sub-001_task-ContinuousVideoGamePlay_run-01_eeg.set
            bids_path = BIDSPath(subject=subject_id,task="ContinuousVideoGamePlay",run="0"+str(run),
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
            raw.set_eeg_reference('average', projection=True)
            mne.set_eeg_reference(raw, ref_channels='average', projection=True)
            
        #Oddball and Gamble dataset
        if(run == 1):

            annotations = raw.annotations 

            # Find all events with "STATUS" in the description
            # Ignore them if they are the last or first event
            status_events = [
                annotation for annotation in annotations
                if "STATUS" in annotation['description'] and annotation != annotations[0] and annotation != annotations[-1]
            ]

            status_times = [annotation['onset'] for annotation in status_events]

            #TODO TODO USE ODDBALL DONE AND ODDBALL START EVENTS, RATHER THAN THIS...

            #This affects subjects 002, 003 and 005
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
                # Subject 006 played 3 sessions of gambling. Maybe just merge them? Or do run 1-3?
                continue
            

            # Split the raw data
            raw_oddball = raw.copy().crop(tmin=raw.times[0], tmax=status_times[0])
            raw_gambling = raw.copy().crop(tmin=status_times[0], tmax=raw.times[-1])

            # For subject 004, the tasks are switched
            if(subject_id == "004"):
                raw_oddball, raw_gambling = raw_gambling, raw_oddball

            #Referencing to average and reconstructing CPz, this is done in mne-bids-pipeline instead
            #raw_oddball, _ = mne.set_eeg_reference(raw_oddball, ref_channels='average')
            #raw_oddball = mne.add_reference_channels(raw_oddball, ref_channels=['CPz'])

            # Writing the oddball task data
            oddball_bids_path = BIDSPath(subject=subject_id, task="Oddball", run="01",
                                        datatype='eeg', suffix='eeg', root=bids_output_root)
            write_raw_bids(raw_oddball, oddball_bids_path, overwrite=True, allow_preload=True, format="EEGLAB")

            # Writing the gambling task data
            gambling_bids_path = BIDSPath(subject=subject_id, task="Gambling", run="01",
                                        datatype='eeg', suffix='eeg', root=bids_output_root)
            write_raw_bids(raw_gambling, gambling_bids_path, overwrite=True, allow_preload=True, format="EEGLAB")


            update_json(oddball_bids_path, "defined in _electrodes.tsv", "AFz", "CPz")
            update_json(gambling_bids_path, "defined in _electrodes.tsv", "AFz", "CPz")
            #TODO events.json also lacks info from original dataset
            #As well as dataset_description.json, participants.json and participants.tsv
            #This is not really a priority. The json files are not directly used for anything in the pipeline (?)
        
        #Axon dataset
        else:
            if(subject_id in ["002","003","005","006","008"]):
                print("Skipping subject " + subject_id)
                continue



            #TODO do this a bit cleaner?
            #Also, dont hardcode the 0.5 second parameter
            new_annotations = [[onset, duration, description] for onset, duration, description in zip(raw.annotations.onset, raw.annotations.duration, raw.annotations.description)]
            annotations = raw.annotations 

            #For Axon, remove superfluous collect star and crash wall events (those occuring within a 500ms window)
            #As it is done in the paper
            for i in range(1,len(annotations)-1):
                for ignores in ["COLLECT_STAR", "PLAYER_CRASH_WALL"]:
                    if(annotations[i]['description'] == ignores):
                        j = 1
                        while(i+j < len(annotations) and (annotations[i+j]['onset'] - annotations[i]['onset'] < 0.5)):
                            #print(i+j, annotations[i+j]['onset'])
                            if(annotations[i+j]['description'] == ignores):
                                new_annotations[i+j][2] = "IGNORE " + ignores
                            j += 1

            
            new_annotations = mne.Annotations(
                onset=[new_annotations[0] for new_annotations in new_annotations],
                duration=[new_annotations[1] for new_annotations in new_annotations],
                description=[new_annotations[2] for new_annotations in new_annotations]
            )

            raw.set_annotations(new_annotations)

            raw.set_eeg_reference('average', projection=True)
            mne.set_eeg_reference(raw, ref_channels='average', projection=True)

            # Don't have to change anything for now
            axon_bids_path = BIDSPath(subject=subject_id, task="Axon", run="01",
                                        datatype='eeg', suffix='eeg', root=bids_output_root)
            write_raw_bids(raw, axon_bids_path, overwrite=True, allow_preload=True, format="EEGLAB")
