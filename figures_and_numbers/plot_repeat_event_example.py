"""Show the Pz signal around repeated wall events in subject 001."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1]


def main():
    path = next((ROOT/'v1.0.0/sub-001/eeg').glob('*run-02_eeg.fdt'))
    channels = path.with_name(path.name.replace('_eeg.fdt', '_channels.tsv'))
    events = ROOT/'prepared_final_bids/sub-001/eeg'/path.name.replace('_eeg.fdt', '_events.tsv')
    with channels.open() as stream:
        names = [r['name'] for r in csv.DictReader(stream, delimiter='\t')]
    with events.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    sf = json.loads(path.with_suffix('.json').read_text())['SamplingFrequency']
    start = 1039.890
    wall = [r for r in rows if r['source_trial_type'] == 'PLAYER_CRASH_WALL'
            and start <= float(r['onset']) < start+.5]
    offsets = np.array([float(r['onset'])-start for r in wall])
    samples = np.memmap(path, mode='r', dtype='<f4').reshape(-1, len(names))
    # Original FDT scalp values are in microvolts, referenced to acquisition CPz.
    signal = sosfiltfilt(butter(4, [.1, 40], btype='bandpass', fs=sf, output='sos'),
                        np.asarray(samples[:, names.index('Pz')], dtype=float))
    left, right = round((start-.3)*sf), round((start+1)*sf)
    values = signal[left:right].copy()
    values -= signal[round((start-.2)*sf):round(start*sf)].mean()
    times = np.arange(left, right)/sf-start
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.plot(times, values, color='#484848', lw=1)
    ax.axvspan(0, .5, color='#d8e7ee', alpha=.6, zorder=-2)
    ax.text(.25, .94, '500 ms eligibility window', ha='center', va='top',
            transform=ax.get_xaxis_transform(), fontsize=9)
    colors = ['#167694']+['#bd623f']*3
    for offset, color in zip(offsets, colors):
        ax.axvline(offset, color=color, lw=.85, alpha=.85)
    ax.set(ylabel='Pz (µV)', ylim=(-16, 17), xlim=(-.3, 1),
           xlabel='Time from first corrected wall event (s)',
           xticks=[-.2, 0, .2, .4, .6, .8, 1],
           title='Subject 001 · run 02 · first wall event at 1039.890 s')
    ax.legend(handles=[Line2D([], [], color='#167694', label='Retained wall event'),
                       Line2D([], [], color='#bd623f', label='Excluded repeat')],
              loc='lower left', ncol=2, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    for extension in ('svg', 'png'):
        fig.savefig(ROOT/f'figures/repeated_wall_example.{extension}', dpi=160)
    plt.close(fig)



if __name__ == '__main__':
    main()
