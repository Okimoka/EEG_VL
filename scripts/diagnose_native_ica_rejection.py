"""Replay only AutoReject on subjects 005/014 to explain saved ICA selection.

No ICA fitting, exclusion application or native-file edits. Reconstructed
retained samples and rejection masks must match the completed cohort.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
os.environ.setdefault('MPLBACKEND','Agg')
import json
from pathlib import Path
import autoreject
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
COHORT=ROOT/'artifacts/ica_native_cohort'
OUT=COHORT/'rejection_diagnostics'

def main():
    OUT.mkdir(exist_ok=True)
    ledger=pd.read_csv(ROOT/'reference/ica_native_cohort/window_diagnostics.tsv',sep='\t',dtype={'subject':str,'run':str})
    results=[]
    for sub in ('005','014'):
        prefix=COHORT/f'sub-{sub}/eeg/sub-{sub}_task-ContinuousVideoGamePlay'
        parts=[]
        for run in ('01','02'):
            raw=mne.io.read_raw_fif(str(prefix)+f'_run-{run}_proc-filt_raw.fif',preload=True,verbose='error')
            raw.filter(1.,100.,n_jobs=1,verbose='error')
            events=mne.make_fixed_length_events(raw,id=3000,duration=4.,overlap=0.,stop=raw.times[-1]-4.)
            part=mne.Epochs(raw,events,event_id={'rest':3000},tmin=0.,tmax=3.996,
                baseline=None,proj=False,preload=True,reject=None,verbose='error')
            parts.append(part)
        epochs=mne.concatenate_epochs(parts,on_mismatch='warn',verbose='error')
        del parts,raw,part
        epochs.set_eeg_reference('average',projection=True,verbose='error')
        epochs.apply_proj(verbose='error')
        saved=mne.read_epochs(str(prefix)+'_proc-icafit_epo.fif',preload=False,verbose='error')
        assert len(epochs)==len(saved.drop_log)
        picks=mne.pick_types(epochs.info,eeg=True,exclude='bads')
        names=[epochs.ch_names[i] for i in picks]
        data=epochs.get_data(picks=picks)
        maximum_difference=0.
        for start in range(0,len(saved),32):
            selected=saved.selection[start:start+32]
            delta=np.abs(data[selected]-saved[start:start+32].get_data(picks=names,verbose='error'))
            maximum_difference=max(maximum_difference,float(delta.max()))
        assert maximum_difference<1e-9,maximum_difference
        detector=autoreject.AutoReject(n_interpolate=[4,8,16],random_state=2026,n_jobs=4,verbose=False)
        detector.fit(epochs)
        reject=detector.get_reject_log(epochs)
        original=np.ones(len(epochs),dtype=bool);original[saved.selection]=False
        assert np.array_equal(reject.bad_epochs,original),f'{sub}: replay disagrees'
        labels=reject.labels[:,picks]
        votes=(labels>0).sum(axis=1)
        ptp=np.ptp(data,axis=-1)*1e6
        fixed=ptp.max(axis=1)>500
        frame=ledger[ledger.subject.eq(sub)].copy().reset_index(drop=True)
        assert frame.retained.to_numpy().tolist()==(~original).tolist()
        frame['bad_channel_votes']=votes
        frame['actual_max_ptp_uv']=ptp.max(axis=1)
        frame['fixed_500uv_reject']=fixed
        frame.to_csv(OUT/f'sub-{sub}_windows.tsv',sep='\t',index=False)
        thresholds=np.array([detector.threshes_[ch]*1e6 for ch in names])
        channel_rows=pd.DataFrame(dict(channel=names,threshold_uv=thresholds,
            fraction_flagged_all=(labels>0).mean(axis=0),
            fraction_flagged_rejected=(labels[original]>0).mean(axis=0),
            median_ptp_uv=np.median(ptp,axis=0),p95_ptp_uv=np.quantile(ptp,.95,axis=0)))
        channel_rows.to_csv(OUT/f'sub-{sub}_channels.tsv',sep='\t',index=False)
        np.savez_compressed(OUT/f'sub-{sub}_reject_log.npz',labels=labels,
            bad_epochs=original,ch_names=np.array(names),thresholds_uv=thresholds)
        groups=[]
        for task,part in frame.groupby('inferred_task'):
            groups.append(dict(task=task,candidates=len(part),autoreject_retained=int(part.retained.sum()),
                              fixed_500uv_retained=int((~part.fixed_500uv_reject).sum()),
                              median_bad_channel_votes=float(part.bad_channel_votes.median()),
                              median_max_ptp_uv=float(part.actual_max_ptp_uv.median())))
        record=dict(subject=sub,identical_rejection_mask=True,max_retained_sample_difference_v=maximum_difference,
            good_eeg_channels=len(names),consensus=float(detector.consensus_['eeg']),
            n_interpolate=int(detector.n_interpolate_['eeg']),
            minimum_bad_channel_votes_in_rejected=int(votes[original].min()),
            median_bad_channel_votes_in_rejected=float(np.median(votes[original])),
            autoreject_rejected=int(original.sum()),fixed_500uv_rejected=int(fixed.sum()),
            rejected_by_both=int((fixed&original).sum()),fixed_only_rejected=int((fixed&~original).sum()),
            threshold_min_uv=float(thresholds.min()),threshold_max_uv=float(thresholds.max()),
            median_max_ptp_rejected_uv=float(np.median(ptp.max(axis=1)[original])),
            median_max_ptp_retained_uv=float(np.median(ptp.max(axis=1)[~original])),
            top_flagged_channels=channel_rows.nlargest(8,'fraction_flagged_rejected').to_dict(orient='records'),
            task_coverage=groups)
        fig,axes=plt.subplots(2,1,figsize=(12,6),sharex=True,gridspec_kw={'height_ratios':[2,1]})
        axes[0].imshow((labels>0).T,aspect='auto',origin='lower',interpolation='nearest',cmap='Greys',vmin=0,vmax=1)
        axes[0].set_yticks(range(len(names)),names,fontsize=5)
        axes[0].set_ylabel('Good EEG channels: black = over learned threshold')
        x=np.arange(len(epochs)); axes[1].scatter(x,votes,c=np.where(original,'#b44a43','#487da4'),s=5)
        axes[1].axhline(detector.consensus_['eeg']*len(names),color='black',ls=':',label='Consensus cutoff')
        axes[1].set(xlabel='Four-second fitting window, ordered by run/time',ylabel='Channels over threshold')
        axes[1].legend(fontsize=8)
        fig.suptitle(f'Subject {sub}: exact AutoReject replay (no ICA refit)\nRed windows rejected; blue retained; no interpolation applied',fontsize=11)
        fig.tight_layout();fig.savefig(OUT/f'sub-{sub}_rejection.png',dpi=160);plt.close(fig)
        # Median-in-time examples avoid selecting only extreme artifacts.
        task='oddball' if sub=='005' else 'gambling'
        chosen=[]
        for rejected in (False,True):
            candidates=np.flatnonzero(frame.inferred_task.eq(task).to_numpy() & (original==rejected))
            chosen.append(int(candidates[len(candidates)//2]))
        fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True)
        example_rows=[]
        for axis,index in zip(axes,chosen):
            y=data[index]*1e6;centered=y-np.median(y,axis=1,keepdims=True)
            axis.plot(epochs.times,centered.T,color='gray',alpha=.3,lw=.5)
            peak=int(ptp[index].argmax());axis.plot(epochs.times,centered[peak],color='#a34035',lw=1,label=names[peak])
            axis.set(title=f'{"Rejected" if original[index] else "Retained"} {task}, run {frame.iloc[index].run}, {frame.iloc[index].start_s:.0f} s; {votes[index]} channels over threshold',ylabel='EEG (µV; individual medians removed)')
            axis.legend(fontsize=8)
            example_rows.append(dict(candidate_index=index,run=frame.iloc[index].run,start_s=float(frame.iloc[index].start_s),rejected=bool(original[index]),votes=int(votes[index]),peak_channel=names[peak],max_ptp_uv=float(ptp[index].max())))
        axes[-1].set_xlabel('Time in window (s)')
        fig.suptitle(f'Subject {sub}: median-in-time retained and rejected {task} examples\nAll good EEG channels; 1–100 Hz, average reference; no component removal',fontsize=11)
        fig.tight_layout();fig.savefig(OUT/f'sub-{sub}_examples.png',dpi=160);plt.close(fig)
        record['examples']=example_rows
        results.append(record)
        (OUT/f'sub-{sub}_summary.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record),flush=True)
        del data,epochs,detector,saved
    (OUT/'summary.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':main()
