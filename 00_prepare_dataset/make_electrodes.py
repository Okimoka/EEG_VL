"""Regenerate electrodes.tsv from the authors' BESA template using MNE."""
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen
import mne
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
URL = 'https://raw.githubusercontent.com/sccn/dipfit/master/standard_BESA/standard-10-5-cap385.elp'

# MNE reads BESA .elp files; remove the first line containing the electrode count.
with TemporaryDirectory() as folder:
    elp = Path(folder) / 'besa.elp'
    elp.write_text('\n'.join(urlopen(URL).read().decode().splitlines()[1:]))
    montage = mne.channels.read_custom_montage(elp, head_size=0.085)

# Match the original dataset's 0.01 mm rounding, keeping template landmarks precise.
positions = montage.get_positions()
montage = mne.channels.make_dig_montage(
    ch_pos={name: np.round(xyz, 5) for name, xyz in positions['ch_pos'].items()},
    nasion=positions['nasion'], lpa=positions['lpa'], rpa=positions['rpa'],
    coord_frame='unknown',
)
# MNE uses the nasion and two ear landmarks to orient and centre the coordinates.
positions = mne.channels.transform_to_head(montage).get_positions()['ch_pos']

# Keep the dataset's channel names and order; all 34 recordings share this layout.
source = HERE.parent / 'v1.0.0/sub-001/eeg/sub-001_task-ContinuousVideoGamePlay_run-01_electrodes.tsv'
table = pd.read_csv(source, sep='\t')
for axis, column in enumerate(('x', 'y', 'z')):
    table[column] = [format(positions[name][axis], '.12g') for name in table['name']]
table.to_csv(HERE / 'electrodes.tsv', sep='\t', index=False, lineterminator='\r\n')
