"""Two selected blink examples; EEG and auxiliary channels have separate scales."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1]


def main():
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharey=True)
    for ax, sub, peak, title in zip(axes, ('001', '014'), (800.298, 344.988),
                                   ('Blink visible in VEOG', 'Blink strongest in HEOG')):
        folder = ROOT/f'v1.0.0/sub-{sub}/eeg'
        names = pd.read_csv(next(folder.glob('*run-01_channels.tsv')), sep='\t').name.tolist()
        raw = np.memmap(next(folder.glob('*run-01_eeg.fdt')), mode='r', dtype='<f4').reshape(-1, 65)
        data = np.array(raw[:, [names.index(n) for n in ('Fp1','Pz','VEOG','HEOG')]], dtype=float).T * 1e-6
        data = sosfiltfilt(butter(4, [.5, 15], btype='bandpass', fs=500, output='sos'), data)
        traces = np.array([data[0]-data[1], data[2], data[3]])
        center = round(peak*500)
        segment = traces[:, center-750:center+750].copy()
        segment -= np.median(segment, axis=1, keepdims=True)
        correlations = np.corrcoef(segment)[0]
        scale = np.percentile(abs(traces[1:, max(0,center-30000):center+30000]), 99)
        times = np.arange(center-750, center+750)/500
        for row, color in enumerate(('#484848', '#2676ae', '#b3693f')):
            sign = -1 if row and correlations[row] < 0 else 1
            ax.plot(times, sign*segment[row]/(100e-6 if row == 0 else scale)+(2-row)*3, color=color, lw=1)
        ax.axvline(peak, color='0.75', ls=':', lw=.8, zorder=-1)
        ax.set(title=f'Subject {sub}\n{title}', xlabel='Time in run 01 (s)',
               xlim=(times[0], times[-1]), ylim=(-.65, 9))
        ax.set_yticks([0,3,6], ['HEOG','VEOG','Fp1−Pz'])
        ax.tick_params(labelleft=True)
        ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    for extension in ('svg', 'png'):
        fig.savefig(ROOT/f'figures/eog_examples.{extension}', dpi=160)


if __name__ == '__main__':
    main()
