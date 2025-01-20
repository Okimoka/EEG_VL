

This repository uses the mne bids pipeline to preprocess and analyze an EEG dataset.
The dataset consists of an Oddball experiment, a gambling Experiment, and an "Axon experiment" (video game) for each of the 17 subjects.


[For collaborators]
The mne bids pipeline allows a workflow where you only need to specify parameters for describing the full process
(the process description is in Pipeline/config.py), and it will automatically apply all of the steps in the correct order, then generate reports.

Some steps could not be covered by mne-bids-pipeline, which I have put into a separate file (01_split_data.py).
The main thing it does is split the gambling experiment and the oddball experiment into two files, so they can be analyzed separately.
But it also does some necessary preprocessing to get the data into a format that mne-bids-pipeline can read.
Afaik the paper applies the same preprocessing on all three experiments.


To get mne-bids-pipeline to work
- install mne_bids_pipeline
- extract dataset ds003517.zip in Pipeline/data such that Pipeline/data contains all the subject folders
- run 01_split_data.py. This should generate files in data_split
- run the mne-bids-pipeline with the command "mne_bids_pipeline config.py --task Oddball" (replace oddball with gambling or axon if desired)
- This took around 2h on my machine for all subjects. You can see my console output in text_output.txt. It should generate reports inside Pipeline/data_split/derivatives
- You can also see my generated results in the .htmls and .tsvs in the repo (e.g. in derivatives/mne-bids-pipeline)


Some things to work on (sorted by priority):
- Perform statistical analysis as described in the "Statistical analyses and classiﬁcation" section of the paper.
Ideally we can identify with statistical significance which game events correspond to which exemplar events.

- Adjust preprocessing params. For example, I think the ICA parameters are not ideal. In the paper they describe

"Thirteen participants had two components for eye movements and four had one
component. In all cases except one, the blink or movement-related components were the ﬁrst ones
(the exception had the 1st and 3rd components)"

Currently, we have much stronger variance in the amount of EOG components. I'm also getting this warning from mne-bids-pipeline:
"No EOG-related ICs detected, this is highly suspicious. A manual check is suggested. You may wish to lower "ica_eog_threshold"."

- I realized too late that the splitting of gambling and oddball data can be done much cleaner.
This will also allow us to take more subjects into consideration which I had simply just skipped for now. 

