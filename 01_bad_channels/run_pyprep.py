"""Run all PyPREP checks with seed 2026 and 1,000 RANSAC subsets."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json
import os
import subprocess
import sys

# One numerical thread per worker; four recordings run in parallel below.
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
import mne
import numpy as np
from pyprep import NoisyChannels

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'logs/pyprep'


def review_windows(extra, channels):
    # Export the failed 1 s / 5 s intervals so the review tool can jump to them.
    windows = []
    for criterion, metric, seconds, threshold in [
        ('bad_by_correlation', 'max_correlations', 1, .4),
        ('bad_by_dropout', 'dropouts', 1, 0),
        ('bad_by_ransac', 'ransac_correlations', 5, .75),
    ]:
        values = np.asarray(extra[criterion][metric])
        # Put every detector's array in channel-by-window order.
        if criterion == 'bad_by_ransac':
            values = values.T
        failed = values > threshold if criterion == 'bad_by_dropout' else values < threshold
        for channel, window in zip(*np.where(failed)):
            windows.append(dict(channel=channels[channel], criterion=criterion,
                start=int(window)*seconds, stop=(int(window)+1)*seconds,
                value=float(values[channel, window]), threshold=threshold))
    return windows


def run(path):
    subject, run = path.name.split('_')[0][4:], path.name.split('run-')[1][:2]
    key = f'sub-{subject}_run-{run}'
    # Read the pipeline's filtered recording; only scalp EEG enters PyPREP.
    raw = mne.io.read_raw_fif(path, preload=True, verbose='error')
    noisy = NoisyChannels(raw.copy().pick('eeg'), do_detrend=True, random_state=2026,
                         ransac=False, correlation=True, reject_by_annotation=None)
    # Run the standard detectors, then RANSAC with our larger sample count.
    noisy.find_all_bads(ransac=False, channel_wise=False)
    # find_all_bads does not expose n_samples, so run RANSAC separately.
    noisy.random_state = np.random.RandomState(2026)
    noisy.find_bad_by_ransac(n_samples=1000, channel_wise=False)
    # Save detector arrays, channel flags and intervals for manual inspection.
    channels = list(noisy.ch_names_original)
    extra = noisy._extra_info
    np.savez_compressed(OUT/f'{key}_diagnostics.npz', **{
        f'{criterion}__{name}': np.asarray(value)
        for criterion, metrics in extra.items() for name, value in metrics.items()})
    result = dict(subject=subject, run=run, recording=key, eeg_ch_names=channels,
        input=str(path.relative_to(ROOT)), settings=dict(random_state=2026, n_samples=1000),
        bads=noisy.get_bads(as_dict=True, verbose=False), all_bads=noisy.get_bads(verbose=False),
        windows=review_windows(extra, channels))
    (OUT/f'{key}.json').write_text(json.dumps(result, indent=2))
    print(key, result['all_bads'], flush=True)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted((ROOT/'logs/filtered').glob('sub-*/eeg/*_proc-filt_raw.fif'))
    with ProcessPoolExecutor(4) as pool:
        list(pool.map(run, files))
    # Write the shared bad-channel flags and event corrections into final BIDS sidecars.
    subprocess.run([sys.executable, str(ROOT/'00_prepare_dataset/prepare_bids.py'), '--final'], check=True)
