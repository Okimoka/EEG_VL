"""Compare subject 007's pilot and final ICA maps saved as small TSV tables."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]


def main():
    tables = [pd.read_csv(ROOT/f'logs/plot_data/sub007_{name}_maps.tsv', sep='\t').set_index('channel')
              for name in ('ica_pilot', 'ica_final_full')]
    common = tables[0].index.intersection(tables[1].index)
    a, b = [table.loc[common].filter(like='ICA').to_numpy() for table in tables]
    correlations = np.corrcoef(a[:, :10].T, b.T)[:10, 10:]
    _, matched = linear_sum_assignment(-abs(correlations))
    fig, axes = plt.subplots(2, 10, figsize=(17.2, 4.75))
    for row, (table, order) in enumerate(zip(tables, (range(10), matched))):
        info = mne.create_info(list(table.index), 250, 'eeg')
        positions = dict(zip(table.index, table[['x', 'y', 'z']].to_numpy()))
        info.set_montage(mne.channels.make_dig_montage(ch_pos=positions, coord_frame='head'))
        for ax, component in zip(axes[row], order):
            weights = table[f'ICA{component:03d}'].to_numpy()
            im, _ = mne.viz.plot_topomap(weights/max(abs(weights)), info, axes=ax, show=False,
                                         cmap='RdBu_r', vlim=(-1, 1), contours=4, res=128)
            ax.set_title(f'ICA{component:03d}')
    fig.subplots_adjust(left=.012, right=.988, bottom=.235, top=.923, wspace=.11, hspace=.35)
    bar = fig.colorbar(im, cax=fig.add_axes([.35, .12, .3, .025]), orientation='horizontal')
    bar.set_label('Normalized component weight')
    for extension in ('svg', 'png'):
        fig.savefig(ROOT/f'figures/sub007_ica_first10_comparison.{extension}', dpi=180)
    print('Final component order:', matched.tolist())


if __name__ == '__main__':
    main()
