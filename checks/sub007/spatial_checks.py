"""Independent descriptive checks of spatial similarity and actual repair burden."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/checks/sub007'
OUT.mkdir(parents=True,exist_ok=True)
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json
import mne,numpy as np,pandas as pd
mne.set_log_level('ERROR')
prefix=ROOT/'artifacts/postica/sub-007/eeg/sub-007_task-ContinuousVideoGamePlay'
epochs=mne.read_epochs(str(prefix)+'_proc-ica_epo.fif',preload=True,verbose='error')
clean=mne.read_epochs(str(prefix)+'_proc-clean_epo.fif',preload=True,verbose='error')
d=np.load(OUT/'diagnostics.npz');names=d['names'].tolist();labels=d['labels'];bad=d['bad_epochs']
x=epochs.get_data(picks=names);xyz=np.array([epochs.info['chs'][epochs.ch_names.index(n)]['loc'][:3] for n in names])
records=[]
for name in ['AFz','F4','F5','FC5','TP9','FC6']:
    index=names.index(name);near=np.argsort(np.linalg.norm(xyz-xyz[index],axis=1))[1:6]
    a=x[:,index].copy();b=x[:,near].mean(axis=1);a-=a.mean(axis=1,keepdims=True);b-=b.mean(axis=1,keepdims=True)
    corr=(a*b).sum(axis=1)/np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
    flagged=labels[:,index]>0
    records.append(dict(channel=name,neighbours=[names[n] for n in near],median_correlation=float(np.median(corr)),median_correlation_flagged=float(np.median(corr[flagged])) if flagged.any() else None,median_correlation_not_flagged=float(np.median(corr[~flagged])),flagged_epochs=int(flagged.sum())))
actual=clean.get_data(picks=names);clean_ptp=np.ptp(actual,axis=-1)*1e6
retained_positions=np.flatnonzero(~bad)
assert np.array_equal(epochs.selection[retained_positions],clean.selection)
rows=pd.read_csv(OUT/'epochs.tsv',sep='\t')
largest=[]
for index in np.argsort(d['ptp_uv'][~bad].max(axis=1))[-5:][::-1]:
    original=retained_positions[index];c=int(np.argmax(d['ptp_uv'][original]));largest.append(dict(trial_id=rows.iloc[original].trial_id,peak_channel_before=names[c],max_ptp_before_uv=float(d['ptp_uv'][original].max()),max_ptp_after_uv=float(clean_ptp[index].max()),peak_channel_repaired=bool(labels[original,c]==2)))
summary=dict(spatial_similarity=records,median_max_ptp_after_native_repairs_uv=float(np.median(clean_ptp.max(axis=1))),largest_original_retained=largest,note='Spatial similarity is descriptive; shared artifacts can also correlate. PTP reduction alone does not establish preservation of neural EEG.')
(OUT/'spatial_checks.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
