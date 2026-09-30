"""Summarize newly fitted native models/windows without importing pilot history."""
from collections import Counter
from pathlib import Path
import re
import mne
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/ica_native_cohort'

def main():
    rows=[];windows=[]
    log=(OUT/'native_pipeline.log').read_text()
    timings={s:float(t) for t,s in re.findall(r'Fitting ICA took ([\d.]+)s\.\n[^\n]*sub-(\d+) Fit',log)}
    for number in range(1,18):
        sub=f'{number:03d}';prefix=OUT/f'sub-{sub}/eeg/sub-{sub}_task-ContinuousVideoGamePlay'
        model=mne.preprocessing.read_ica(str(prefix)+'_proc-ica_ica.fif',verbose='error')
        epochs=mne.read_epochs(str(prefix)+'_proc-icafit_epo.fif',preload=False,verbose='error')
        table=pd.read_csv(str(prefix)+'_proc-ica_components.tsv',sep='\t')
        good=[epochs.ch_names[i] for i in mne.pick_types(epochs.info,eeg=True,exclude='bads')]
        assert good==model.ch_names and model.n_components_==len(good)-1
        assert epochs.info['sfreq']==250 and len(epochs.times)==1000 and epochs.baseline is None
        assert np.allclose([epochs.info['highpass'],epochs.info['lowpass']],[1,100])
        assert model.fit_params['extended'] and model.method=='infomax'
        labels=Counter()
        for index in model.exclude:
            text=table.loc[table.component.eq(index),'status_description'].iloc[0]
            match=re.match(r'Auto-detected (.+) \(MNE-ICALabel\)',str(text))
            if not match:raise ValueError(f'{sub}/{index}: missing original ICLabel description')
            labels[match[1].replace(' ','_')]+=1
        rows.append(dict(subject=sub,candidate_windows=len(epochs.drop_log),retained_windows=len(epochs),
            rejected_windows=len(epochs.drop_log)-len(epochs),rejected_percent=100*(1-len(epochs)/len(epochs.drop_log)),
            retained_minutes=len(epochs)*4/60,components=model.n_components_,proposed_exclusions=len(model.exclude),
            current_status_bad=int(table.status.eq('bad').sum()),fit_seconds=timings.get(sub,float('nan')),
            **{'proposed_'+key:labels[key] for key in ('eye_blink','muscle_artifact','heart_beat','line_noise','channel_noise')}))
        index=0
        for run in ('01','02'):
            raw=mne.io.read_raw_fif(str(prefix)+f'_run-{run}_proc-filt_raw.fif',preload=False,verbose='error')
            grid=mne.make_fixed_length_events(raw,duration=4.,overlap=0.,stop=raw.times[-1]-4.)
            for event in grid:
                windows.append(dict(subject=sub,candidate=index,run=run,start_s=(event[0]-raw.first_samp)/250,
                    retained=index in set(epochs.selection)))
                index+=1
        assert index==len(epochs.drop_log)
    pd.DataFrame(rows).to_csv(OUT/'results_summary.tsv',sep='\t',index=False)
    pd.DataFrame(windows).to_csv(OUT/'window_diagnostics.tsv',sep='\t',index=False)
    print(pd.DataFrame(rows).sum(numeric_only=True)[['candidate_windows','retained_windows','components','proposed_exclusions']])

if __name__=='__main__':main()
