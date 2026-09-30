"""Read-only replay of subject 007 post-ICA AutoReject. Local outputs only."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/checks/sub007'
OUT.mkdir(parents=True,exist_ok=True)
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ['MPLBACKEND']='Agg'
os.environ['JOBLIB_TEMP_FOLDER']=str(OUT/'joblib_tmp')
os.environ['MPLCONFIGDIR']=str(OUT/'mpl_cache')
import hashlib,json,time
import autoreject,joblib,mne
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    started=time.time();mne.set_log_level('ERROR')
    prefix=ROOT/'artifacts/postica/sub-007/eeg/sub-007_task-ContinuousVideoGamePlay'
    paths=[Path(str(prefix)+s) for s in ('_proc-ica_epo.fif','_proc-clean_epo.fif','_proc-ica_components.tsv','_proc-ica_ica.fif')]
    paths += [ROOT/'config_postica.py',ROOT/'artifacts/postica/trial_ledger.tsv']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    (OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    epochs=mne.read_epochs(paths[0],preload=True,verbose='error')
    clean=mne.read_epochs(paths[1],preload=False,verbose='error')
    picks=mne.pick_types(epochs.info,eeg=True,exclude='bads')
    names=[epochs.ch_names[i] for i in picks]
    values=epochs.get_data(picks=picks)*1e6
    ptp=np.ptp(values,axis=-1)
    original=~np.isin(epochs.selection,clean.selection)
    ledger=pd.read_csv(paths[-1],sep='\t',dtype={'subject':str,'run':str})
    rows=ledger[ledger.subject.eq('007') & ledger.boundary_retained].copy().reset_index(drop=True)
    assert rows.native_candidate.tolist()==epochs.selection.tolist()
    assert np.array_equal(rows.final_retained.to_numpy(),~original)
    rows['task']=np.where(rows.condition.str.startswith('ODDBALL'),'oddball',np.where(rows.condition.str.startswith('GAMBLING'),'gambling','gameplay'))
    # Choose examples before replay or any alternate rejection/ERP inspection.
    selected=[]
    for task in ('oddball','gambling','gameplay'):
        fig,axes=plt.subplots(3,2,figsize=(12,8),sharex=True)
        for col,bad in enumerate((False,True)):
            candidates=np.flatnonzero(rows.task.eq(task).to_numpy() & (original==bad))
            for row,q in enumerate((.25,.5,.75)):
                index=int(candidates[int((len(candidates)-1)*q)])
                peak=int(ptp[index].argmax());y=values[index]-np.median(values[index],axis=1,keepdims=True)
                ax=axes[row,col];ax.plot(epochs.times,y.T,color='gray',alpha=.3,lw=.5)
                ax.plot(epochs.times,y[peak],color='#a34035',lw=1,label=names[peak]);ax.axvline(0,color='k',lw=.5,ls=':')
                ax.set(title=f'{"Rejected" if bad else "Retained"} · {rows.iloc[index].corrected_onset_s:.3f} s · {ptp[index].max():.0f} µV PTP',ylabel='EEG (µV)');ax.title.set_fontsize(9);ax.legend(fontsize=7)
                selected.append(dict(task=task,index=index,rejected=bool(bad),chronological_quantile=q,trial_id=rows.iloc[index].trial_id,max_ptp_uv=float(ptp[index].max()),peak_channel=names[peak]))
        for ax in axes[-1]:ax.set_xlabel('Time relative to event (s)')
        fig.suptitle(f'007 {task}: predetermined trace sample, condition labels hidden\nAfter ICA, before local interpolation; individual channel medians removed',fontsize=11)
        fig.tight_layout();fig.savefig(OUT/f'{task}_initial_examples.png',dpi=145,bbox_inches='tight');plt.close(fig)
    (OUT/'initial_examples.json').write_text(json.dumps(selected,indent=2)+'\n')
    print('Initial traces saved; beginning exact AutoReject replay.',flush=True)
    detector=autoreject.AutoReject(n_interpolate=np.array([4,8,16]),random_state=2026,n_jobs=4,verbose=False)
    detector.fit(epochs)
    reject=detector.get_reject_log(epochs)
    match=bool(np.array_equal(reject.bad_epochs,original))
    labels=reject.labels[:,picks]
    votes=(labels>0).sum(axis=1)
    thresholds=np.array([detector.threshes_[name]*1e6 for name in names])
    assert match,'Replayed mask does not match native output; stop sensitivity analysis.'
    joblib.dump(detector,OUT/'autoreject_replay.joblib')
    np.savez_compressed(OUT/'diagnostics.npz',labels=labels,ptp_uv=ptp,thresholds_uv=thresholds,names=names,bad_epochs=original,votes=votes,selection=epochs.selection)
    rows['bad_channel_votes']=votes;rows['max_ptp_uv']=ptp.max(axis=1)
    rows.to_csv(OUT/'epochs.tsv',sep='\t',index=False)
    channels=pd.DataFrame(dict(channel=names,threshold_uv=thresholds,flagged_percent=100*(labels>0).mean(axis=0),repaired_retained_percent=100*(labels[~original]==2).mean(axis=0),median_ptp_uv=np.median(ptp,axis=0)))
    channels.to_csv(OUT/'channels.tsv',sep='\t',index=False)
    records=[]
    for condition,group in rows.groupby('condition',sort=False):
        ids=group.index.to_numpy();bad=original[ids]
        records.append(dict(condition=condition,total=len(ids),retained=int((~bad).sum()),median_votes=float(np.median(votes[ids])),median_ptp_uv=float(np.median(ptp[ids].max(axis=1)))))
    fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True,gridspec_kw={'height_ratios':[2,1]})
    axes[0].imshow((labels>0).T,aspect='auto',origin='lower',interpolation='nearest',cmap='Greys')
    axes[0].set_yticks(np.arange(len(names)),names,fontsize=5);axes[0].set_ylabel('Channels over learned thresholds')
    axes[1].scatter(np.arange(len(rows)),votes,c=np.where(original,'#b75538','#147d92'),s=5)
    axes[1].axhline(detector.consensus_['eeg']*len(names),color='k',ls=':');axes[1].set(xlabel='Epoch in time order',ylabel='Channel votes')
    fig.suptitle('007: exact post-ICA rejection replay; red rejected, blue retained',fontsize=11);fig.tight_layout();fig.savefig(OUT/'rejection_votes.png',dpi=160);plt.close(fig)
    for p in paths:assert digest(p)==hashes[str(p.relative_to(ROOT))]
    result=dict(exact_native_mask=match,epochs=len(rows),retained=int((~original).sum()),good_eeg=len(names),consensus=float(detector.consensus_['eeg']),n_interpolate=int(detector.n_interpolate_['eeg']),minimum_votes_rejected=int(votes[original].min()),median_votes_rejected=float(np.median(votes[original])),median_votes_retained=float(np.median(votes[~original])),threshold_range_uv=[float(thresholds.min()),float(thresholds.max())],conditions=records,wall_seconds=time.time()-started,input_hashes_unchanged=True)
    (OUT/'baseline_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
