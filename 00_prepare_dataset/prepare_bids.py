"""
Prepare ds003517 to create a dataset that can be read by mne-bids-pipeline

use --final after PyPREP to add bads and event corrections.
"""
from pathlib import Path
import json
import shutil
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'v1.0.0'
TEMPLATES = ROOT / '00_prepare_dataset'
EVENTS = {'ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING WIN', 'GAMBLING LOSS',
          'PLAYER_CRASH_WALL', 'COLLECT_STAR', 'MISSILE_HIT_ENEMY',
          'PLAYER_CRASH_ENEMY', 'COLLECT_AMMO', 'SHOOT_BUTTON'}


def correct_events(table):
    last = {}
    result = []
    for number, row in enumerate(table.to_dict('records'), 1):
        label = row['trial_type']
        sample = round(float(row['onset']) * 500)
        # 500 ms rule: compare each wall/star event with the last retained event of its type.
        keep = True
        if label in ('PLAYER_CRASH_WALL', 'COLLECT_STAR'):
            keep = sample - last.get(label, -1000000) >= 250
            if keep:
                last[label] = sample
        # Visual delay: move analysis events 40 ms later
        shift = 0.04 if label in EVENTS else 0
        # Keep the original timing and labels in extra columns.
        row.update(source_event_row=str(number), source_onset=row['onset'],
                   source_sample=row['sample'], source_trial_type=label,
                   onset_shift_seconds=f'{shift:.3f}', keep_same_type_500ms=str(keep).lower(),
                   exclusion_reason='n/a' if keep else 'same-type event <500 ms after last retained')
        row['onset'] = f"{float(row['onset']) + shift:.10f}"
        if row['sample'] != 'n/a':
            row['sample'] = str(int(row['sample']) + round(shift * 500))
        # Rename repeated events so epoch selection skips them
        if not keep:
            row['trial_type'] = 'EXCLUDED_REPEAT__' + label
        result.append(row)
    return pd.DataFrame(result)


def prepare(final=False):
    destination = ROOT / ('prepared_final_bids' if final else 'prepared_bids')
    # Copy editable sidecars; link the unchanged signals.
    # Replace existing links first so metadata edits never reach the originals.
    for source in SOURCE.rglob('*'):
        if source.is_file() and source.suffix != '.lock':
            target = destination / source.relative_to(SOURCE)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.unlink(missing_ok=True)
            if source.suffix in ('.tsv', '.json'):
                shutil.copyfile(source, target)
            else:
                target.symlink_to(source.resolve())
    for folder in sorted(destination.glob('sub-*/eeg')):
        subject = folder.parent.name
        # Both runs need the same bad-channel set for the shared ICA fit.
        if final:
            bads = [json.loads((ROOT/f'logs/pyprep/{subject}_run-{r}.json').read_text())['all_bads'] for r in ('01', '02')]
            union = set(bads[0] + bads[1])
        for run, path in enumerate(sorted(folder.glob('*_channels.tsv'))):
            table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
            # Fix the original n/a channel types and units: 63 EEG and two EOG.
            table['type'] = ['EOG' if name in ('VEOG', 'HEOG') else 'EEG' for name in table.name]
            table['units'] = 'uV'
            if final:
                table['status'] = ['bad' if name in union else 'good' for name in table.name]
                table['status_description'] = ['PyPREP union of both runs' if name in union else 'n/a' for name in table.name]
                table['pyprep_bad_in_this_run'] = [str(name in bads[run]).lower() for name in table.name]
                path.with_suffix('.json').write_text(json.dumps({'pyprep_bad_in_this_run': {'Description': 'Per-run PyPREP flag before the subject union.'}}, indent=2))
            table.to_csv(path, sep='\t', index=False)
        # Reuse the corrected tables (identical for all runs).
        # electrodes.tsv was generated using make_electrodes.py
        for name in ('electrodes.tsv', 'coordsystem.json'):
            for path in folder.glob('*_' + name):
                shutil.copyfile(TEMPLATES/name, path)
        # Make the recording JSON agree with the corrected channel types.
        for path in folder.glob('*_eeg.json'):
            info = json.loads(path.read_text())
            info.update(EEGChannelCount=63, EOGChannelCount=2, MiscChannelCount=0,
                        DigitizedLandmarks=False, DigitizedHeadPoints=False)
            path.write_text(json.dumps(info, indent=2))
        # Apply the event corrections only in the final prepared dataset.
        if final:
            for path in folder.glob('*_events.tsv'):
                table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
                correct_events(table).to_csv(path, sep='\t', index=False)
    print(destination.name)


if __name__ == '__main__':
    prepare(final='--final' in sys.argv)
