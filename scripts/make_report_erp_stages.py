"""Plot saved participant ERPs without changing preprocessing or native results."""
import os
import argparse
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import matplotlib.pyplot as plt
import mne
import numpy as np
from _plotting import mean_ci

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'artifacts/postica'
OUTPUT = ROOT / 'report/figures'
AUDIT = ROOT / 'artifacts/report_restructure'
CONDITIONS = ('ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING WIN', 'GAMBLING LOSS')

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    global SOURCE
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reference',action='store_true');args=parser.parse_args()
    if args.reference: SOURCE=ROOT/'reference/erp'
    mne.set_log_level('ERROR')
    OUTPUT.mkdir(parents=True,exist_ok=True)
    AUDIT.mkdir(parents=True,exist_ok=True)
    data = {stage: defaultdict(list) for stage in ('Before_ICA', 'Final_cleaned')}
    participants = defaultdict(list)
    hashes = {}
    for number in range(1, 18):
        sub = f'{number:03d}'
        pairs = {}
        for stage in data:
            path = SOURCE / (f'evokeds/sub-{sub}_{stage}_ave.fif' + ('.gz' if args.reference else ''))
            hashes[str(path.relative_to(ROOT))] = digest(path)
            pairs[stage] = {e.comment: e for e in mne.read_evokeds(path, verbose='error')}
        for label in CONDITIONS:
            if label not in pairs['Before_ICA']:
                assert label not in pairs['Final_cleaned']
                continue
            before, final = (pairs[stage][label] for stage in data)
            assert before.nave == final.nave
            assert before.ch_names == final.ch_names
            assert np.array_equal(before.times, final.times)
            assert np.allclose(before.baseline, [-.2, 0]) and np.allclose(final.baseline, [-.2, 0])
            assert not before.info['bads'] and not final.info['bads']
            for stage in data:
                data[stage][label].append(pairs[stage][label])
            participants[label].append(sub)
    assert [len(participants[c]) for c in CONDITIONS] == [17, 17, 14, 14]
    # Check that saved before-ICA values match the original matched-ERP plot data.
    for label in CONDITIONS:
        original = np.load(SOURCE / f'grand_{label.replace(" ", "_")}.npz')
        channel = 'Pz' if label.startswith('ODDBALL') else 'Cz'
        observed = np.array([e.get_data(picks=[channel])[0] for e in data['Before_ICA'][label]]) * 1e6
        assert list(original['subjects']) == participants[label]
        assert np.allclose(observed, original['Before_ICA'], rtol=1e-6, atol=1e-5)
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.4), sharex=True, sharey='col')
    colors = {'Standard': '#2878ad', 'Rare': '#c45136', 'Win': '#2878ad', 'Loss': '#c45136'}
    for row, (stage, title) in enumerate((('Before_ICA', 'Before ICA'), ('Final_cleaned', 'ICA + AutoReject repairs'))):
        for col, (labels, channel, task) in enumerate(((CONDITIONS[:2], 'Pz', 'Oddball'), (CONDITIONS[2:], 'Cz', 'Gambling'))):
            ax = axes[row, col]
            evokeds = {label.split()[-1].title(): data[stage][label] for label in labels}
            mne.viz.plot_compare_evokeds(evokeds, picks=channel, axes=ax, ci=mean_ci,
                colors={label: colors[label] for label in evokeds}, show_sensors=False,
                show=False, vlines=[0], truncate_yaxis=False, truncate_xaxis=False, legend='upper left',
                title=f'{title}\n{task} · {channel} (N={len(participants[labels[0]])})')
            ax.set_xlim(-.2, .8)
            ax.set_xticks([-.2, 0, .2, .4, .6, .8])
            ax.tick_params(labelbottom=True)
            ax.set_xlabel('Time after corrected visual onset (s)')
    for col in range(2):
        values = []
        labels = CONDITIONS[:2] if col == 0 else CONDITIONS[2:]
        channel = 'Pz' if col == 0 else 'Cz'
        for stage in data:
            for label in labels:
                values.extend(mean_ci(np.array([e.get_data(picks=[channel])[0] for e in data[stage][label]]) * 1e6).ravel())
        low, high = min(values), max(values)
        margin = .12 * (high - low)
        axes[0, col].set_ylim(low - margin, high + margin)
    fig.suptitle('Grand-average condition ERPs · 0.1–40 Hz', fontsize=14)
    fig.text(.5, .015, 'Same final retained trials in both rows; equal subject weights.\nShading: pointwise 95% t confidence intervals. Matched reference, global display interpolation and baseline.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .075, 1, .96))
    for extension in ('svg', 'pdf', 'png'):
        fig.savefig(OUTPUT / f'erp_cleaning_stages.{extension}', dpi=180)
    plt.close(fig)
    assert all(digest(ROOT / name) == expected for name, expected in hashes.items())
    (AUDIT / 'erp_figure_validation.json').write_text(json.dumps(dict(
        passed=True, input_sha256=hashes, participants=participants,
        matched_retained_trial_counts=True, before_ica_matches_original_plot=True,
        input_files_unchanged=True, baseline=[-.2, 0], band_hz=[.1, 40],
        caveat='Bottom row includes both ICA removal and final local AutoReject repairs.'), indent=2)+'\n')
    print('Created four-panel ERP figure; saved inputs and matched trial counts verified.')

if __name__ == '__main__':
    main()
