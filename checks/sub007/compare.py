"""Two predeclared fixed-rejection comparisons; no ICA refitting or main writes."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/checks/sub007'
OUT.mkdir(parents=True,exist_ok=True)
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(OUT/'mpl_cache')
import hashlib,json
import mne,numpy as np,pandas as pd
import matplotlib.pyplot as plt

def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    baseline=json.loads((OUT/'baseline_summary.json').read_text());assert baseline['exact_native_mask']
    mne.set_log_level('ERROR')
    d=np.load(OUT/'diagnostics.npz');rows=pd.read_csv(OUT/'epochs.tsv',sep='\t',dtype={'subject':str,'run':str})
    names=d['names'].tolist();ptp=d['ptp_uv'];labels=d['labels'];bad=d['bad_epochs'];votes=d['votes'];thresholds=d['thresholds_uv']
    prefix=ROOT/'artifacts/postica/sub-007/eeg/sub-007_task-ContinuousVideoGamePlay'
    epochs=mne.read_epochs(str(prefix)+'_proc-ica_epo.fif',preload=True,verbose='error')
    clean=mne.read_epochs(str(prefix)+'_proc-clean_epo.fif',preload=True,verbose='error')
    assert np.array_equal(epochs.selection,d['selection'])
    values=epochs.get_data(picks=names)*1e6
    keep={'Native AutoReject':~bad,'Fixed 75 µV':ptp.max(axis=1)<=75,'Fixed 100 µV':ptp.max(axis=1)<=100}
    records=[];metrics=[];chosen=[]
    for variant,mask in keep.items():
        added=mask&bad;lost=~mask&~bad
        for condition,g in rows.groupby('condition',sort=False):
            ids=g.index.to_numpy();records.append(dict(variant=variant,condition=condition,available=len(ids),retained=int(mask[ids].sum()),recovered_vs_native=int(added[ids].sum()),lost_vs_native=int(lost[ids].sum())))
        metrics.append(dict(variant=variant,retained=int(mask.sum()),rejected=int((~mask).sum()),recovered_vs_native=int(added.sum()),lost_vs_native=int(lost.sum()),median_max_ptp_before_repairs_uv=float(np.median(ptp[mask].max(axis=1))),median_votes_before_repairs=float(np.median(votes[mask])),median_votes_recovered=float(np.median(votes[added])) if added.any() else None,local_interpolation=variant=='Native AutoReject'))
        if variant=='Native AutoReject':continue
        tag='75' if '75' in variant else '100'
        candidates=np.flatnonzero(added & rows.task.eq('oddball').to_numpy())
        fig,axes=plt.subplots(3,1,figsize=(11,7),sharex=True)
        for ax,q in zip(axes,(.25,.5,.75)):
            index=int(candidates[int((len(candidates)-1)*q)])
            y=values[index]-np.median(values[index],axis=1,keepdims=True);peak=int(ptp[index].argmax())
            ax.plot(epochs.times,y.T,color='gray',alpha=.3,lw=.5);ax.plot(epochs.times,y[peak],color='#a34035',lw=1,label=names[peak]);ax.axvline(0,color='k',ls=':',lw=.5)
            ax.set(ylabel='EEG (µV)',title=f'Onset {rows.iloc[index].corrected_onset_s:.3f} s · maximum {ptp[index].max():.1f} µV · {votes[index]} channels above native thresholds');ax.title.set_fontsize(9);ax.legend(fontsize=8)
            chosen.append(dict(variant=variant,index=index,trial_id=rows.iloc[index].trial_id,maximum_uv=float(ptp[index].max()),votes=int(votes[index]),peak_channel=names[peak]))
        axes[-1].set_xlabel('Time relative to event (s)');fig.suptitle(f'007: oddball epochs recovered by {variant}; chronological quartiles\nNo local repairs; condition labels hidden; individual medians removed',fontsize=11)
        fig.tight_layout();fig.savefig(OUT/f'recovered_{tag}_examples.png',dpi=150,bbox_inches='tight');plt.close(fig)
    pd.DataFrame(records).to_csv(OUT/'comparison_counts.tsv',sep='\t',index=False)
    (OUT/'comparison_summary.json').write_text(json.dumps(metrics,indent=2)+'\n')
    (OUT/'recovered_examples.json').write_text(json.dumps(chosen,indent=2)+'\n')
    for variant,mask in keep.items():rows[variant]=mask
    rows.to_csv(OUT/'comparison_trial_ledger.tsv',sep='\t',index=False)
    # Inspect targeted frequently repaired channels in the same preset oddball examples.
    initial=json.loads((OUT/'initial_examples.json').read_text());initial=[e for e in initial if e['task']=='oddball']
    good_info=mne.pick_info(epochs.info,mne.pick_channels(epochs.ch_names,names,ordered=True))
    locs=np.array([c['loc'][:3] for c in good_info['chs']]);centers=['AFz','F4','F5']
    neighbour_records=[]
    for center in centers:
        c=names.index(center);near=np.argsort(np.linalg.norm(locs-locs[c],axis=1))[1:6]
        fig,axes=plt.subplots(3,2,figsize=(12,7.5),sharex=True)
        for j,example in enumerate(initial):
            index=example['index'];ax=axes[j%3,j//3];selected=[c]+near.tolist();y=values[index,selected].copy();y-=np.median(y,axis=1,keepdims=True)
            for k,n in enumerate(selected):ax.plot(epochs.times,y[k],lw=1 if k==0 else .65,color='#a34035' if k==0 else 'gray',alpha=1 if k==0 else .6,label=names[n] if k==0 else None)
            ax.set(ylabel='EEG (µV)',title=f'{"Rejected" if bad[index] else "Retained"} · {rows.iloc[index].corrected_onset_s:.3f} s\n{center}: {ptp[index,c]:.1f} µV PTP; threshold {thresholds[c]:.1f} µV');ax.title.set_fontsize(9);ax.legend(fontsize=8)
            neighbour_records.append(dict(channel=center,trial_id=rows.iloc[index].trial_id,ptp_uv=float(ptp[index,c]),threshold_uv=float(thresholds[c]),neighbours=[names[n] for n in near]))
        for ax in axes[-1]:ax.set_xlabel('Time relative to event (s)')
        fig.suptitle(f'007: {center} (red) and five nearest good EEG neighbours (gray)\nSame preset oddball sample, before local repairs; each channel median removed',fontsize=11);fig.tight_layout();fig.savefig(OUT/f'{center}_neighbours.png',dpi=150,bbox_inches='tight');plt.close(fig)
    (OUT/'neighbour_examples.json').write_text(json.dumps(neighbour_records,indent=2)+'\n')
    # Matched individual epochs: show how native repair changes these examples.
    fig,axes=plt.subplots(3,2,figsize=(12,7.5),sharex=True)
    original_baselined=epochs.copy().apply_baseline((-.2,0),verbose='error')
    for j,example in enumerate(initial[:3]):
        index=example['index'];ci=int(np.flatnonzero(clean.selection==epochs.selection[index])[0])
        for k,ch in enumerate(['AFz','Pz']):
            ax=axes[j,k];before=original_baselined[index].get_data(picks=[ch])[0,0]*1e6;after=clean[ci].get_data(picks=[ch])[0,0]*1e6
            ax.plot(epochs.times,before,label='ICA only',color='#147d92',lw=1);ax.plot(epochs.times,after,label='ICA + local AR',color='#b75538',lw=1)
            ax.set(ylabel='EEG (µV)',title=f'{ch} · retained oddball · onset {rows.iloc[index].corrected_onset_s:.3f} s');ax.title.set_fontsize(9);ax.legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('Time relative to event (s)')
    fig.suptitle('007: native repair effect on identical retained epochs\nBoth baseline corrected; read with the stored average-reference projector',fontsize=11);fig.tight_layout();fig.savefig(OUT/'matched_repair_examples.png',dpi=150,bbox_inches='tight');plt.close(fig)
    channels=pd.read_csv(OUT/'channels.tsv',sep='\t');channels['ptp_to_threshold_median']=channels.median_ptp_uv/channels.threshold_uv
    channels.sort_values('flagged_percent',ascending=False).to_csv(OUT/'channels_ranked.tsv',sep='\t',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    x=np.arange(len(names));axes[0].plot(x,thresholds,label='Learned threshold',lw=1);axes[0].plot(x,np.median(ptp,axis=0),label='Median observed PTP',lw=1);axes[0].set_xticks(x,names,rotation=90,fontsize=6);axes[0].set_ylabel('µV');axes[0].legend(fontsize=8)
    for task in ['oddball','gambling','gameplay']:
        idx=rows.task.eq(task).to_numpy();axes[1].hist(votes[idx],bins=np.arange(0,len(names)+2)-.5,density=True,alpha=.4,label=task)
    axes[1].axvline(baseline['consensus']*len(names),color='k',ls=':',label='Rejection boundary');axes[1].set(xlabel='Channels above their thresholds',ylabel='Fraction per bin');axes[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'thresholds_and_votes.png',dpi=160,bbox_inches='tight');plt.close(fig)
    hashes=json.loads((OUT/'input_hashes.json').read_text())
    for path,expected in hashes.items():assert digest(ROOT/path)==expected,path
    (OUT/'validation.json').write_text(json.dumps(dict(passed=True,main_inputs_unchanged=True,no_ica_refit=True,no_erp_or_classification_optimization=True,variants=list(keep),input_hashes=hashes),indent=2)+'\n')
    print(json.dumps(metrics,indent=2));print(channels.sort_values('flagged_percent',ascending=False).head(10).to_string(index=False))

if __name__=='__main__':main()
