"""Compare regular/Picard Infomax for subject 001 on identical training samples."""
from pathlib import Path
import json
import mne
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
ROOT=Path(__file__).resolve().parents[1]

def main():
    models=[];epochs=[]
    for name in ('ica_native_cohort','ica_picard_pilot'):
        prefix=ROOT/f'artifacts/{name}/sub-001/eeg/sub-001_task-ContinuousVideoGamePlay'
        models.append(mne.preprocessing.read_ica(str(prefix)+'_proc-ica_ica.fif',verbose='error'))
        epochs.append(mne.read_epochs(str(prefix)+'_proc-icafit_epo.fif',preload=False,verbose='error'))
    a,b=epochs
    assert a.ch_names==b.ch_names and np.array_equal(a.selection,b.selection)
    assert np.array_equal(a.events,b.events) and np.array_equal(a.times,b.times)
    for start in range(0,len(a),32):assert np.array_equal(a[start:start+32].get_data(),b[start:start+32].get_data())
    assert models[0].ch_names==models[1].ch_names
    maps=[m.get_components() for m in models]
    norm=[]
    for x in maps:
        x=x-x.mean(axis=0);norm.append(x/np.linalg.norm(x,axis=0))
    correlation=norm[0].T@norm[1]
    left,right=linear_sum_assignment(-np.abs(correlation))
    rows=[dict(infomax=int(i),picard=int(j),absolute_map_correlation=float(abs(correlation[i,j])),
        infomax_automatic_exclusion=i in models[0].exclude,picard_automatic_exclusion=j in models[1].exclude) for i,j in zip(left,right)]
    out=ROOT/'artifacts/ica_picard_pilot';pd.DataFrame(rows).to_csv(out/'matched_components.tsv',sep='\t',index=False)
    print(json.dumps(dict(identical_fitting_samples=True,matched_maps_over_0p998=sum(r['absolute_map_correlation']>.998 for r in rows),
        automatic_exclusions_agree=all(r['infomax_automatic_exclusion']==r['picard_automatic_exclusion'] for r in rows)),indent=2))
    print('Compare the native "Fitting ICA took ..." lines in both logs; runtimes depend on hardware.')

if __name__=='__main__':main()
