"""
Prepare cleaned EEG trials and random gameplay controls.
Partly assisted by LLM
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import autoreject
import joblib
import mne
import numpy as np
import pandas as pd

POSTICA = Path('../artifacts/postica')
WORK = Path('../artifacts/analysis_final')
mne.set_log_level('ERROR')


def features(epochs):
    # Convert epochs to 63 EEG channels at 250 Hz, in microvolts.
    epochs.pick('eeg').interpolate_bads()
    epochs.resample(250, method='polyphase').crop(0, .796)
    return epochs.get_data().astype('float32') * 1e6


def prepare(subject, event_rows):
    folder = WORK / f'sub-{subject}'
    folder.mkdir(parents=True, exist_ok=True)
    prefix = POSTICA / f'sub-{subject}/eeg/sub-{subject}_task-ContinuousVideoGamePlay'
    native = mne.read_epochs(f'{prefix}_proc-clean_epo.fif', preload=True)
    # Keep the exemplary trials and the three gameplay event types.
    wanted = ['ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING LOSS', 'GAMBLING WIN',
              'MISSILE_HIT_ENEMY', 'PLAYER_CRASH_ENEMY', 'PLAYER_CRASH_WALL']
    native = native[native.metadata.event_name.isin(wanted).to_numpy()]

    # Recover the AutoReject settings used to clean the gameplay controls.
    ica = mne.read_epochs(f'{prefix}_proc-ica_epo.fif', preload=True)
    ar = autoreject.AutoReject(n_interpolate=[4, 8, 16], random_state=2026, n_jobs=2, verbose=False).fit(ica)
    joblib.dump(ar, folder / 'autoreject.joblib')
    ar.get_reject_log(ica).save(folder / 'source_reject_log.npz', overwrite=True)
    np.save(folder / 'source_epoch_indices.npy', ica.selection)

    # Randomly phased 2-second grids inside gameplay spans, split at 30-second gaps.
    candidates, raws = [], {}
    gameplay = event_rows[~event_rows.condition.str.startswith(('ODDBALL', 'GAMBLING'))]
    for run, rows in gameplay.groupby('run'):
        raw = mne.io.read_raw_fif(f'{prefix}_run-{run}_proc-clean_raw.fif')
        raws[run] = raw
        samples = np.sort(rows.corrected_sample.unique())
        blocks = np.split(samples, np.flatnonzero(np.diff(samples) > 15000) + 1)
        rng = np.random.default_rng(2026 + int(subject) * 100 + int(run))
        for block in blocks:
            start, stop = max(int(block[0]), 100), min(int(block[-1]), raw.n_times - 401)
            if stop - start >= 500:
                offset = int(rng.integers(0, 1000))
                candidates.extend((run, sample) for sample in range(start + offset, stop, 1000))
    rng = np.random.default_rng(2026 + int(subject))
    chosen = rng.choice(len(candidates), 400, replace=False)
    selected = pd.DataFrame(sorted(candidates[i] for i in chosen), columns=['run', 'sample'])
    selected.to_csv(folder / 'control_candidates.tsv', sep='\t', index=False)
    parts = []
    for run, rows in selected.groupby('run'):
        events = np.column_stack([rows['sample'], np.zeros(len(rows), int), np.full(len(rows), 999)])
        parts.append(mne.Epochs(raws[run], events, event_id={'RANDOM_GAMEPLAY': 999},
            tmin=-.2, tmax=.8, baseline=None, metadata=rows, preload=True))

    # Repair or reject control epochs before applying the baseline correction.
    controls = mne.concatenate_epochs(parts)
    clean, log = ar.transform(controls, return_log=True)
    log.save(folder / 'control_reject_log.npz', overwrite=True)
    controls.metadata.assign(retained=~log.bad_epochs).to_csv(folder / 'controls.tsv', sep='\t', index=False)
    clean.apply_baseline((-.2, 0)).save(folder / 'controls-epo.fif', overwrite=True)
    # Read back the float32 FIF, as for the native event epochs.
    clean = mne.read_epochs(folder / 'controls-epo.fif', preload=True)
    condition = np.r_[native.metadata.event_name.to_numpy(), np.repeat('RANDOM_GAMEPLAY', len(clean))]
    eeg = np.concatenate([features(native), features(clean)])
    np.savez_compressed(WORK / f'sub-{subject}.npz', eeg=eeg, condition=condition.astype(str))
    print(subject, flush=True)


if __name__ == '__main__':
    events = pd.read_csv(POSTICA / 'eligible_events.tsv', sep='\t', dtype={'subject': str, 'run': str})
    for subject in [f'{i:03d}' for i in range(1, 18)]:
        prepare(subject, events[events.subject.eq(subject)])
