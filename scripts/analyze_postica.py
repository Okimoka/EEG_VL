"""Validate native post-ICA outputs and plot matched EEG/ERP comparisons.

Only reads native results. No fitting, rejection, reclassification or pipeline
patching. Display copies interpolate global bads and baseline identically.
The ICA comparison uses the final retained IDs but no per-epoch AR repairs.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from _plotting import mean_ci

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'artifacts/postica'
ERP = ['ODDBALL STANDARD', 'ODDBALL RARE', 'GAMBLING WIN', 'GAMBLING LOSS']

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')

def display(epochs, indices, label, eeg):
    evoked = epochs[indices].average(picks=eeg)
    if evoked.info['bads']:
        evoked.interpolate_bads(reset_bads=True, verbose='error')
    evoked.apply_baseline((-.2, 0), verbose='error')
    evoked.comment = label
    return evoked

def plot_erps(grand, subjects):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4))
    colors = {'Before ICA': '#147d92', 'After ICA': '#b75538'}
    for ax, label in zip(axes.flat, ERP):
        channel = 'Pz' if label.startswith('ODDBALL') else 'Cz'
        data = {stage: grand[stage][label] for stage in colors}
        mne.viz.plot_compare_evokeds(data, picks=channel, axes=ax, ci=mean_ci,
            colors=colors, show_sensors=False, show=False, vlines=[0],
            truncate_yaxis=False, title=f'{label.title()} · {channel} (N={len(subjects[label])})')
        ax.set_xlabel('Time after corrected visual onset (s)')
        arrays = {stage.replace(' ', '_'): np.array([e.copy().pick(channel).data[0] for e in values])*1e6 for stage, values in data.items()}
        np.savez_compressed(OUT/f'grand_{label.replace(" ", "_")}.npz',
            times=data['Before ICA'][0].times, subjects=subjects[label], channel=channel, **arrays)
    fig.suptitle('Matched ERPs before and after ICA · 0.1–40 Hz', fontsize=14)
    fig.text(.5, .02, 'Identical final retained trials, reference, display interpolation and baseline. No per-epoch repairs in either curve.\nShading: pointwise 95% confidence interval across equally weighted participants.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .06, 1, .95))
    for ext in ('png', 'pdf', 'svg'): fig.savefig(OUT/f'erp_ica_comparison.{ext}', dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, labels, channel in zip(axes, (ERP[:2], ERP[2:]), ('Pz', 'Cz')):
        mne.viz.plot_compare_evokeds({label: grand['Final cleaned'][label] for label in labels},
            picks=channel, axes=ax, ci=mean_ci, show_sensors=False, show=False,
            vlines=[0], truncate_yaxis=False, title=f'Final cleaned ERPs · {channel}')
        ax.set_xlabel('Time after corrected visual onset (s)')
    fig.suptitle('ICA + post-ICA AutoReject · 0.1–40 Hz', fontsize=13)
    fig.tight_layout()
    for ext in ('png', 'pdf'): fig.savefig(OUT/f'final_cleaned_erps.{ext}', dpi=170)
    plt.close(fig)

def main():
    mne.set_log_level('ERROR')
    preflight = json.loads((OUT/'preflight.json').read_text())
    ledger = pd.read_csv(OUT/'eligible_events.tsv', sep='\t', dtype={'subject': str, 'run': str})
    policy = json.loads((ROOT/'artifacts/review/channel_policy.json').read_text())
    log = (OUT/'native_pipeline.log').read_text()
    assert 'pipeline exit status: 0' in log
    assert 'Fitting ICA took' not in log and 'Fitting ICA to data' not in log
    ar_matches = re.findall(r"sub-(\d+) autoreject marked (\d+) epochs as bad .*?\{'eeg': (?:np\.int\d+\()?([0-9]+)", log)
    ar_parameters = {s: int(n) for s, _, n in ar_matches}
    ar_removed = {s: int(n) for s, n, _ in ar_matches}
    assert len(ar_parameters) == 17
    protected = {str(p.relative_to(ROOT)): digest(p) for p in OUT.glob('sub-*/eeg/*') if p.is_file() and not p.name.endswith('.lock')}
    grand = {stage: defaultdict(list) for stage in ('Before ICA', 'After ICA', 'Final cleaned')}
    participants = defaultdict(list)
    counts, summaries, validations, channel_repairs, all_trials = [], [], [], [], []
    waveform_dir = OUT/'evokeds'; waveform_dir.mkdir(exist_ok=True)
    for number in range(1, 18):
        sub = f'{number:03d}'; prefix = OUT/f'sub-{sub}/eeg/sub-{sub}_task-ContinuousVideoGamePlay'
        path = lambda suffix: Path(str(prefix)+suffix)
        before = mne.read_epochs(path('_epo.fif'), preload=False, verbose='error')
        after = mne.read_epochs(path('_proc-ica_epo.fif'), preload=False, verbose='error')
        clean = mne.read_epochs(path('_proc-clean_epo.fif'), preload=False, verbose='error')
        events = ledger[ledger.subject.eq(sub)].copy().reset_index(drop=True)
        # MNE omits annotations outside the signal range before epoching.
        # In 012/run-02 the last corrected event is at n_times, just past the file.
        full_event_count = len(events)
        native_keys = {}
        raw_lengths = {}
        for run in ('01', '02'):
            raw = mne.io.read_raw_fif(path(f'_run-{run}_proc-filt_raw.fif'), preload=False, verbose='error')
            raw_lengths[run] = raw.n_times
            ids = {name: i+1 for i, name in enumerate(events.condition.unique())}
            native, ids = mne.events_from_annotations(raw, event_id=ids, verbose='error')
            inverse = {value: key for key, value in ids.items()}
            for sample, _, code in native:
                key = (run, int(sample), inverse[code])
                assert key not in native_keys
                native_keys[key] = len(native_keys)
        events['native_candidate'] = [native_keys.get((row.run, row.corrected_sample, row.condition), -1)
            for row in events.itertuples()]
        unavailable = events[events.native_candidate.lt(0)].copy()
        assert all(row.corrected_sample < 0 or row.corrected_sample >= raw_lengths[row.run]
            for row in unavailable.itertuples())
        unavailable = unavailable.assign(ica_only_max_ptp_uv=np.nan, ica_only_peak_channel='',
            boundary_retained=False, final_retained=False, drop_reason='OUTSIDE_RECORDING',
            changed_channels_after_ar=0)
        events = events[events.native_candidate.ge(0)].sort_values('native_candidate').reset_index(drop=True)
        assert events.native_candidate.tolist() == list(range(len(native_keys)))
        assert len(events) == len(before.drop_log) == len(clean.drop_log)
        assert np.array_equal(before.selection, after.selection)
        assert set(clean.selection) <= set(before.selection)
        assert len(before)-len(clean) == ar_removed[sub]
        assert before.baseline is None and after.baseline is None
        assert np.allclose(clean.baseline, [-.2, 0])
        assert before.info['sfreq'] == after.info['sfreq'] == clean.info['sfreq'] == 500
        assert len(before.times) == 501 and np.allclose(before.times[[0, -1]], [-.2, .8])
        assert np.allclose([before.info['highpass'], before.info['lowpass']], [.1, 40], rtol=0, atol=1e-7)
        assert set(before.info['bads']) == set(clean.info['bads']) == set(policy['subject_union_bads'][sub])
        eeg = [before.ch_names[i] for i in mne.pick_types(before.info, eeg=True, exclude=[])]
        good = [name for name in eeg if name not in before.info['bads']]
        assert len(eeg) == 63
        events['ica_only_max_ptp_uv'] = np.nan
        events['ica_only_peak_channel'] = ''
        for start in range(0, len(after), 32):
            values = after[start:start+32].get_data(picks=good)
            assert np.isfinite(values).all()
            ptp = np.ptp(values, axis=-1)*1e6
            selections = after.selection[start:start+len(values)]
            events.loc[selections, 'ica_only_max_ptp_uv'] = ptp.max(axis=1)
            events.loc[selections, 'ica_only_peak_channel'] = [good[i] for i in ptp.argmax(axis=1)]
        native_names = [dict((value, key) for key, value in before.event_id.items())[code] for code in before.events[:, 2]]
        assert native_names == events.iloc[before.selection].condition.tolist()
        assert before.metadata.event_name.tolist() == native_names
        # Concatenation changes absolute event samples by one fixed offset per run.
        for run in ('01', '02'):
            selected = events.iloc[before.selection].run.eq(run).to_numpy()
            offsets = before.events[selected, 0] - events.iloc[before.selection[selected]].corrected_sample.to_numpy()
            assert len(np.unique(offsets)) == 1
            if run == '01': assert offsets[0] == 0
        model = mne.preprocessing.read_ica(path('_proc-ica_ica.fif'), verbose='error')
        decisions = pd.read_csv(path('_proc-ica_components.tsv'), sep='\t')
        model.exclude = decisions.loc[decisions.status.eq('bad'), 'component'].astype(int).tolist()
        assert model.exclude == preflight['component_exclusions'][sub] and model.ch_names == good
        assert np.allclose([ch['loc'][:3] for ch in model.info['chs']],
                           [before.info['chs'][before.ch_names.index(name)]['loc'][:3] for name in good], atol=1e-8, rtol=0)
        indices = np.unique([0, len(before)//2, len(before)-1])
        probe = before[indices].load_data()
        prediction = model.apply(probe.copy(), verbose='error')
        application_error = float(np.max(np.abs(prediction.get_data(picks=good)-after[indices].get_data(picks=good))))
        assert application_error < 5e-9, application_error
        reference_error = float(np.max(np.abs(probe.get_data(picks=good).mean(axis=1))))
        assert reference_error < 1e-9
        events['boundary_retained'] = events.native_candidate.isin(before.selection)
        events['final_retained'] = events.native_candidate.isin(clean.selection)
        events['drop_reason'] = ['; '.join(reason) for reason in clean.drop_log]
        events['changed_channels_after_ar'] = 0
        repairs = np.zeros(len(good), dtype=int)
        interpolation_count = np.zeros(len(clean), dtype=int)
        # Measure actual repairs relative to ICA-only data with the identical baseline.
        # 1 nV tolerance exceeds FIF rounding error; globally bad EEG is excluded.
        post_positions = np.searchsorted(after.selection, clean.selection)
        for start in range(0, len(clean), 32):
            block = after[post_positions[start:start+32]].load_data().apply_baseline((-.2, 0), verbose='error')
            actual = clean[start:start+32].get_data(picks=good)
            assert np.isfinite(actual).all()
            delta = actual-block.get_data(picks=good)
            # Reading active average-reference projectors distributes an interpolation
            # change over all good channels. Remove that shared shift before counting
            # repaired channels; fewer than half can be repaired (at most 16).
            delta -= np.median(delta, axis=1, keepdims=True)
            changed = np.max(np.abs(delta), axis=-1) > 1e-9
            assert np.all(changed.sum(axis=1) <= ar_parameters[sub])
            repairs += changed.sum(axis=0)
            interpolation_count[start:start+len(changed)] = changed.sum(axis=1)
        events.loc[clean.selection, 'changed_channels_after_ar'] = interpolation_count
        for name, value in zip(good, repairs):
            channel_repairs.append(dict(subject=sub, channel=name, retained_epochs_changed=int(value), retained_epochs=len(clean)))
        saved_evokeds = {stage: [] for stage in grand}
        matched_before = np.searchsorted(before.selection, clean.selection)
        for label, group in events.groupby('condition', sort=False):
            kept = group[group.final_retained]
            assert len(kept), f'{sub}/{label}: no final trials'
            row = dict(subject=sub, condition=label, eligible=len(group)+int(unavailable.condition.eq(label).sum()),
                boundary_retained=int(group.boundary_retained.sum()), final_retained=len(kept),
                postica_rejected=int(group.boundary_retained.sum())-len(kept),
                retained_with_repairs=int((kept.changed_channels_after_ar>0).sum()),
                mean_repaired_channels=float(kept.changed_channels_after_ar.mean()))
            counts.append(row)
            clean_indices = np.flatnonzero(events.iloc[clean.selection].condition.eq(label).to_numpy())
            for stage, epochs, positions in (
                ('Before ICA', before, matched_before[clean_indices]),
                ('After ICA', after, post_positions[clean_indices]),
                ('Final cleaned', clean, clean_indices)):
                evoked = display(epochs, positions, label, eeg)
                saved_evokeds[stage].append(evoked)
                grand[stage][label].append(evoked)
            participants[label].append(sub)
        for stage, evokeds in saved_evokeds.items():
            mne.write_evokeds(waveform_dir/f'sub-{sub}_{stage.replace(" ", "_")}_ave.fif', evokeds, overwrite=True, verbose='error')
        raw_eog_checks = []
        for run in ('01', '02'):
            filt = mne.io.read_raw_fif(path(f'_run-{run}_proc-filt_raw.fif'), preload=False, verbose='error')
            cleaned = mne.io.read_raw_fif(path(f'_run-{run}_proc-clean_raw.fif'), preload=False, verbose='error')
            from mne_bids import BIDSPath, read_raw_bids
            original_aux = read_raw_bids(BIDSPath(root=ROOT/'prepared_native_bids', subject=sub, run=run, task='ContinuousVideoGamePlay', datatype='eeg', suffix='eeg', extension='.set'), verbose='error')
            assert filt.n_times == cleaned.n_times == original_aux.n_times
            for start in range(0, filt.n_times, 100000):
                end = min(start+100000, filt.n_times)
                x = filt.get_data(picks=['VEOG','HEOG'], start=start, stop=end)
                assert np.array_equal(x, cleaned.get_data(picks=['VEOG','HEOG'], start=start, stop=end))
                assert np.allclose(x, original_aux.get_data(picks=['VEOG','HEOG'], start=start, stop=end), rtol=1e-7, atol=1e-12)
            raw_eog_checks.append(run)
        summaries.append(dict(subject=sub, eligible=full_event_count, boundary_retained=len(before), final_retained=len(clean),
            postica_rejected=len(before)-len(clean), rejected_percent=100*(1-len(clean)/len(before)),
            retained_with_repairs=int((interpolation_count>0).sum()),
            repaired_percent=100*float((interpolation_count>0).mean()),
            median_changed_channels=float(np.median(interpolation_count)),
            max_changed_channels=int(interpolation_count.max()), selected_n_interpolate=ar_parameters.get(sub),
            excluded_components=len(model.exclude)))
        validations.append(dict(subject=sub, original_trial_ids_match=True, corrected_event_samples_match=True,
            model_application_max_error_v=application_error, good_channel_reference_error_v=reference_error,
            eog_unchanged_runs=raw_eog_checks, baseline_after_rejection=True, matched_comparison_has_no_local_repairs=True))
        all_trials.append(pd.concat([events, unavailable]).sort_values('candidate'))
        print(f'{sub}: {len(clean)}/{len(before)} retained, {(interpolation_count>0).sum()} repaired epochs', flush=True)
    frame = pd.DataFrame(summaries)
    frame.to_csv(OUT/'results_summary.tsv', sep='\t', index=False)
    pd.DataFrame(counts).to_csv(OUT/'condition_counts.tsv', sep='\t', index=False)
    pd.concat(all_trials).to_csv(OUT/'trial_ledger.tsv', sep='\t', index=False)
    pd.DataFrame(channel_repairs).to_csv(OUT/'channel_repairs.tsv', sep='\t', index=False)
    for stage, conditions in grand.items():
        group_evokeds = []
        for label, values in conditions.items():
            evoked = mne.grand_average(values, interpolate_bads=False)
            evoked.comment = label
            group_evokeds.append(evoked)
        mne.write_evokeds(OUT/f'group_{stage.replace(" ", "_")}_ave.fif', group_evokeds, overwrite=True, verbose='error')
    plot_erps(grand, participants)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].bar(frame.subject, frame.rejected_percent, color='#b75538')
    axes[0].set_ylabel('Analysis epochs rejected (%)')
    axes[1].bar(frame.subject, frame.repaired_percent, color='#147d92')
    axes[1].set(ylabel='Retained epochs repaired (%)', xlabel='Subject')
    fig.suptitle('Post-ICA AutoReject on event-related epochs')
    fig.tight_layout()
    for ext in ('png', 'pdf'): fig.savefig(OUT/f'postica_rejection.{ext}', dpi=160)
    plt.close(fig)
    totals = {key: int(frame[key].sum()) for key in ('eligible', 'boundary_retained', 'final_retained', 'postica_rejected', 'retained_with_repairs', 'excluded_components')}
    totals['rejected_percent'] = 100*totals['postica_rejected']/totals['boundary_retained']
    totals['repaired_percent_of_retained'] = 100*totals['retained_with_repairs']/totals['final_retained']
    summary = dict(totals=totals, subjects=summaries, erp_participants={label: participants[label] for label in ERP}, analysis_participants=dict(participants),
        comparison='Final retained trial IDs, identical global interpolation/baseline; neither curve contains per-epoch AutoReject repairs.',
        final_cleaned='Separate native epochs and final cleaned averages include AutoReject repairs.',
        repair_measure='Channels changed by >1 nV versus baseline-matched ICA-only epochs after removing the channel-common reference shift; globally bad EEG excluded.')
    write_json(OUT/'results_summary.json', summary)
    for name, expected in preflight['protected_native_outputs'].items():
        assert digest(ROOT/name) == expected, f'Changed fitted cohort input: {name}'
    for name, expected in protected.items():
        assert digest(ROOT/name) == expected, f'Analysis changed native output: {name}'
    write_json(OUT/'validation.json', dict(passed=True, checked_utc=datetime.now(timezone.utc).isoformat(), subjects=validations,
        input_ica_cohort_unchanged=True, native_analysis_outputs_unchanged_during_validation=True,
        native_output_sha256=protected, no_refitting_or_reclassification=True,
        script_sha256=digest(Path(__file__)), totals=totals))
    print(json.dumps(totals, indent=2))

if __name__ == '__main__':
    main()
