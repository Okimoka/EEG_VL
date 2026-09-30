"""Stage saved ICA models and decisions for a native post-ICA or broadband pass.

The default uses the exact reviewed models. Fresh fits require --models refit
and --reviewed after manual review of their own component-status TSVs. Component
indices from the reference run are never transferred onto a different fit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import mne
import numpy as np
import pandas as pd
from _prepared_metadata import ANALYSIS_LABELS

ROOT=Path(__file__).resolve().parents[1]

def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch',choices=['postica','ica_broadband'],default='postica')
    parser.add_argument('--models',choices=['reference','refit'],default='reference')
    parser.add_argument('--reviewed',action='store_true',help='Confirm that the freshly fitted component-status tables have been reviewed')
    args=parser.parse_args()
    if args.models=='refit' and not args.reviewed:
        parser.error('Review the newly fitted ICA reports and TSVs, then pass --reviewed.')
    source=ROOT/('reference/ica' if args.models=='reference' else 'artifacts/ica_native_cohort')
    out=ROOT/'artifacts'/args.branch
    policy=json.loads((ROOT/'artifacts/review/channel_policy.json').read_text())
    manifest=json.loads((ROOT/'reference/file_manifest.json').read_text())
    records=[];copies=[];protected={};exclusions={};plans=[]
    # Validate the full cohort before writing any destination models.
    for number in range(1,18):
        sub=f'{number:03d}';prefix=f'sub-{sub}/eeg/sub-{sub}_task-ContinuousVideoGamePlay'
        modelpath=source/(prefix+'_proc-ica_ica.fif')
        decisionspath=source/(prefix+'_proc-ica_components.tsv')
        model=mne.preprocessing.read_ica(modelpath,verbose='error')
        decisions=pd.read_csv(decisionspath,sep='\t')
        assert decisions.component.tolist()==list(range(model.n_components_))
        assert set(decisions.status)<={'good','bad'}
        exclusions[sub]=decisions.loc[decisions.status.eq('bad'),'component'].astype(int).tolist()
        candidate=0
        for run in ('01','02'):
            folder=ROOT/f'prepared_native_bids/sub-{sub}/eeg'
            stem=f'sub-{sub}_task-ContinuousVideoGamePlay_run-{run}'
            channels=pd.read_csv(folder/(stem+'_channels.tsv'),sep='\t',keep_default_na=False)
            good=channels.loc[channels.type.eq('EEG') & channels.status.eq('good'),'name'].tolist()
            assert good==model.ch_names,f'{sub}: prepared channels/reference differ from the model'
            assert sorted(channels.loc[channels.status.eq('bad'),'name'])==policy['subject_union_bads'][sub]
            positions=pd.read_csv(folder/(stem+'_electrodes.tsv'),sep='\t').set_index('name')
            assert np.allclose(positions.loc[good,['x','y','z']], [ch['loc'][:3] for ch in model.info['chs']],atol=1e-8,rtol=0)
            events=pd.read_csv(folder/(stem+'_events.tsv'),sep='\t',keep_default_na=False)
            selected=events[events.trial_type.isin(ANALYSIS_LABELS)]
            assert not selected.onset.duplicated().any()
            for row in selected.itertuples():
                assert row.keep_same_type_500ms in (True,'true')
                assert abs(float(row.onset)-float(row.source_onset)-.04)<1e-8
                source_row=int(row.source_event_row)-1
                records.append(dict(subject=sub,run=run,candidate=candidate,
                    trial_id=f'sub-{sub}_run-{run}_event-{source_row:06d}',source_row_index=source_row,
                    condition=row.trial_type,original_onset_s=float(row.source_onset),
                    corrected_onset_s=float(row.onset),corrected_sample=round(float(row.onset)*500)))
                candidate+=1
        for suffix in ('_proc-ica_ica.fif','_proc-ica_components.tsv'):
            src=source/(prefix+suffix);dst=out/(prefix+suffix);sha=digest(src)
            if args.models=='reference':assert sha==manifest[str(src.relative_to(ROOT))],src
            if dst.exists() and digest(dst)!=sha:
                raise SystemExit(f'Destination contains a different model/decision table: {dst}; use a separate output folder.')
            plans.append((src,dst,sha));protected[str(src.relative_to(ROOT))]=sha
        # Fresh reports retain the native ICLabel table/component plots. Reference
        # replay creates application reports; enormous training HTMLs are omitted.
        if args.models=='refit':
            for suffix in ('_report.html','_report.h5'):
                src=source/(prefix+suffix);dst=out/(prefix+suffix)
                if src.exists() and not dst.exists():plans.append((src,dst,digest(src)))
    for src,dst,sha in plans:
        dst.parent.mkdir(parents=True,exist_ok=True)
        if not dst.exists():shutil.copy2(src,dst)
        copies.append(dict(source=str(src.relative_to(ROOT)),destination=str(dst.relative_to(ROOT)),sha256=sha))
    pd.DataFrame(records).to_csv(out/'eligible_events.tsv',sep='\t',index=False)
    receipt=dict(models=args.models,protected_native_outputs=protected,copied_inputs=copies,
        component_exclusions=exclusions,total_component_exclusions=sum(map(len,exclusions.values())),
        candidate_events=len(records),prepared_events_already_corrected=True,no_ica_refit=True)
    (out/'preflight.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f"Staged {len(exclusions)} models; {receipt['total_component_exclusions']} exclusions; {len(records)} eligible anchors.")

if __name__=='__main__':main()
