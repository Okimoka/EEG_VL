"""
Train EEGNet classifiers and apply them to gameplay events.
Partly assisted by LLM
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from braindecode import EEGClassifier
from braindecode.models import EEGNet
from sklearn.metrics import balanced_accuracy_score
from sklearn.utils.class_weight import compute_class_weight

WORK = Path('../artifacts/analysis_final')
EVENTS = ['MISSILE_HIT_ENEMY', 'PLAYER_CRASH_ENEMY', 'PLAYER_CRASH_WALL']
TASKS = {'oddball': ['ODDBALL STANDARD', 'ODDBALL RARE'], 'gambling': ['GAMBLING LOSS', 'GAMBLING WIN']}
CONTRASTS = ['direct', 'positive_vs_background', 'negative_vs_background']
torch.use_deterministic_algorithms(True)
torch.backends.cudnn.benchmark = False
torch.set_num_threads(2)


def source(data, task, contrast, window):
    # Select either two exemplary categories or one category versus gameplay.
    arrays, labels, people = [], [], []
    for subject, saved in data.items():
        if task == 'gambling' and subject in ['002', '003', '005']:
            continue
        conditions = saved['condition']
        negative, positive = TASKS[task]
        if contrast == 'direct':
            indices = np.flatnonzero(np.isin(conditions, [negative, positive]))
        else:
            positive = positive if contrast == 'positive_vs_background' else negative
            indices = np.r_[np.flatnonzero(conditions == positive), np.flatnonzero(conditions == 'RANDOM_GAMEPLAY')]
        arrays.append(saved['eeg'][indices, :, window])
        # Label 1 is the selected exemplary category; label 0 is the alternative.
        labels.append((conditions[indices] == positive).astype('int64'))
        people.extend([subject] * len(indices))
    return np.concatenate(arrays), np.concatenate(labels), np.array(people)


def fit(x, y, people, data, profile, task, contrast, window, subject):
    # Hold out one participant; calculate scaling and class weights from the others.
    test = people == subject
    mean = x[~test].mean(axis=(0, 2), dtype=np.float64)[None, :, None].astype('float32')
    std = x[~test].std(axis=(0, 2), dtype=np.float64)[None, :, None].astype('float32')
    weights = compute_class_weight('balanced', classes=np.array([0, 1]), y=y[~test])
    target = np.isin(data[subject]['condition'], EVENTS + ['RANDOM_GAMEPLAY'])
    conditions = data[subject]['condition'][target]
    tx = data[subject]['eeg'][target, :, window]
    source_rows, transfer_rows = [], []

    # Train three fresh models with fixed seeds.
    for repeat in range(3):
        seed = (102026 + 200000 * (profile == 'fixed_window')
                + 1000 * (task == 'gambling') + 10000 * CONTRASTS.index(contrast) + int(subject) + repeat)
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        net = EEGClassifier(EEGNet(n_chans=63, n_outputs=2, n_times=x.shape[-1], sfreq=250),
            optimizer=torch.optim.Adam, lr=.001, optimizer__weight_decay=0.,
            criterion__weight=torch.tensor(weights, dtype=torch.float32, device='cuda'),
            max_epochs=20, batch_size=64, train_split=None, iterator_train__drop_last=False, device='cuda', verbose=0)
        net.fit((x[~test] - mean) / std, y[~test])
        # Predict the held-out source trials and gameplay events.
        sp = net.predict_proba((x[test] - mean) / std)[:, 1]
        tp = net.predict_proba((tx - mean) / std)[:, 1]
        torch.save(dict(weights=net.module_.state_dict(), mean=mean, std=std, seed=seed,
            source_probability=sp, target_probability=tp),
            WORK / 'models' / f'{profile}_{task}_{contrast}_{subject}_{repeat}.pt')
        identity = dict(profile=profile, task=task, contrast=contrast, subject=subject, repeat=repeat)
        source_rows.append(dict(**identity, accuracy=balanced_accuracy_score(y[test], sp >= .5)))

        # Count labels assigned at the 0.5 threshold, not average confidence scores.
        controls = tp[conditions == 'RANDOM_GAMEPLAY']
        for event in EVENTS:
            values = tp[conditions == event]
            if min(len(values), len(controls)) >= 20:
                a, b = (values >= .5).mean(), (controls >= .5).mean()
                transfer_rows.append(dict(**identity, event=event, event_fraction=a, control_fraction=b, difference=a-b))
    print(profile, task, contrast, subject, flush=True)
    return source_rows, transfer_rows


if __name__ == '__main__':
    (WORK / 'models').mkdir(parents=True, exist_ok=True)
    data = {f'{i:03d}': dict(np.load(WORK / f'sub-{i:03d}.npz')) for i in range(1, 18)}
    source_rows, transfer_rows = [], []
    for profile in ['whole_epoch', 'fixed_window']:
        for task in TASKS:
            # Sample indices at 250 Hz: 0–800, 300–600 or 200–350 ms.
            window = slice(0, 200) if profile == 'whole_epoch' else (slice(75, 150) if task == 'oddball' else slice(50, 88))
            for contrast in CONTRASTS[:1] if profile == 'whole_epoch' else CONTRASTS:
                x, y, people = source(data, task, contrast, window)
                for subject in np.unique(people):
                    s, t = fit(x, y, people, data, profile, task, contrast, window, subject)
                    source_rows.extend(s)
                    transfer_rows.extend(t)
                pd.DataFrame(source_rows).to_csv(WORK / 'source_runs.tsv', sep='\t', index=False)
                pd.DataFrame(transfer_rows).to_csv(WORK / 'transfer_runs.tsv', sep='\t', index=False)
