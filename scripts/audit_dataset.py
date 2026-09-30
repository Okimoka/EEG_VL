"""Recording/event inventory and unmodified EOG amplitude diagnostics."""
from pathlib import Path
import json
import mne
from mne_bids import BIDSPath,read_raw_bids
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/qc'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    inventory=[];event_counts=[];eog=[]
    for number in range(1,18):
        subject=f'{number:03d}'
        for run in ('01','02'):
            raw=read_raw_bids(BIDSPath(root=ROOT/'prepared_bids',subject=subject,run=run,
                task='ContinuousVideoGamePlay',datatype='eeg',suffix='eeg',extension='.set'),verbose='error')
            events=pd.read_csv(next((ROOT/f'v1.0.0/sub-{subject}/eeg').glob(f'*run-{run}_events.tsv')),sep='\t')
            key=f'sub-{subject}_run-{run}'
            inventory.append(dict(recording=key,subject=subject,run=run,samples=raw.n_times,sfreq=raw.info['sfreq'],
                duration_s=raw.n_times/raw.info['sfreq'],eeg=len(mne.pick_types(raw.info,eeg=True,exclude=[])),
                eog=len(mne.pick_types(raw.info,eog=True,exclude=[])),event_rows=len(events),
                gambling_outcomes=int(events.trial_type.isin(['GAMBLING WIN','GAMBLING LOSS']).sum())))
            for label,count in events.trial_type.value_counts().items():
                event_counts.append(dict(recording=key,subject=subject,run=run,trial_type=label,count=count))
            for ch in ('VEOG','HEOG'):
                x=raw.get_data(picks=[ch])[0]*1e6
                eog.append(dict(recording=key,channel=ch,median_uV=float(np.median(x)),sd_uV=float(x.std()),
                    abs_p99_uV=float(np.percentile(np.abs(x),99)),median_centered_abs_p99_uV=float(np.percentile(np.abs(x-np.median(x)),99)),
                    peak_to_peak_uV=float(np.ptp(x)),max_abs_time_s=float(np.argmax(np.abs(x))/raw.info['sfreq'])))
            print(key,flush=True)
    for name,rows in [('recordings',inventory),('event_counts',event_counts),('eog_metrics',eog)]:
        pd.DataFrame(rows).to_csv(OUT/f'{name}.tsv',sep='\t',index=False)
    summary={ch:{metric:dict(min=float(min(r[metric] for r in eog if r['channel']==ch)),
        median=float(np.median([r[metric] for r in eog if r['channel']==ch])),
        max=float(max(r[metric] for r in eog if r['channel']==ch)))
        for metric in ('median_uV','median_centered_abs_p99_uV','peak_to_peak_uV')} for ch in ('VEOG','HEOG')}
    (OUT/'eog_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
