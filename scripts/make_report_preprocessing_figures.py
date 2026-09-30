"""Render report-only figures from accepted counts and unchanged saved EEG."""
import os
import argparse
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLBACKEND','Agg')
from pathlib import Path
import hashlib,json
import matplotlib.pyplot as plt
import mne,numpy as np,pandas as pd
from _montage import topomap_sphere
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'report/figures';AUDIT=ROOT/'artifacts/report_notes_revision'

def save(fig,name):
    fig.savefig(OUT/f'{name}.svg',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reference',action='store_true');args=parser.parse_args()
    base=ROOT/('reference' if args.reference else 'artifacts')
    mne.set_log_level('ERROR');OUT.mkdir(parents=True,exist_ok=True);AUDIT.mkdir(parents=True,exist_ok=True)
    counts=pd.read_csv(base/'qc/bad_channel_frequency.tsv',sep='\t')
    if args.reference:
        info=mne.io.read_info(base/'plot_data/scalp-info.fif',verbose='error')
    else:
        raw=mne.io.read_raw_fif(ROOT/'artifacts/mne-bids-pipeline/sub-001/eeg/sub-001_task-ContinuousVideoGamePlay_run-01_proc-filt_raw.fif',preload=False,verbose='error')
        info=mne.pick_info(raw.info,mne.pick_types(raw.info,eeg=True,exclude=[]))
    counts=counts.set_index('channel').loc[info.ch_names]
    assert int(counts.all_recordings.sum())==68
    vmax=int(counts.all_recordings.max());fig,axes=plt.subplots(1,3,figsize=(11.5,4.0))
    for ax,column,title in zip(axes,['all_recordings','run_01','run_02'],['All recordings (34)','Run 01: oddball / gambling (17)','Run 02: gameplay (17)']):
        names=[f'{name}\n{int(counts.loc[name,column])}' if counts.loc[name,'all_recordings']>=5 else '' for name in info.ch_names]
        im,_=mne.viz.plot_topomap(counts[column].to_numpy(float),info,axes=ax,names=names,
             contours=0,image_interp='linear',cmap='YlOrRd',vlim=(0,vmax),show=False,sphere=topomap_sphere(),extrapolate='head')
        for text in ax.texts:
            if text.get_text():text.set_bbox(dict(facecolor='white',alpha=.8,edgecolor='none',pad=.8))
        ax.set_title(title,fontsize=11)
    fig.subplots_adjust(left=.02,right=.89,top=.87,bottom=.06,wspace=.2)
    cax=fig.add_axes([.92,.2,.018,.55]);bar=fig.colorbar(im,cax=cax);bar.set_label('Recordings with a flag');bar.set_ticks(range(vmax+1))
    save(fig,'bad_channel_topomap')
    summary=pd.read_csv(base/'ica_native_cohort/results_summary.tsv',sep='\t',dtype={'subject':str})
    x=np.arange(len(summary));median=float(summary.rejected_percent.median())
    fig,ax=plt.subplots(figsize=(10,3.2));ax.bar(x,summary.rejected_percent,color='#597c9c')
    ax.axhline(median,color='black',linestyle=':',linewidth=1,label=f'Median over all subjects ({median:.1f}%)')
    ax.set(ylabel='ICA fitting windows rejected (%)',xlabel='Subject',ylim=(0,55));ax.set_xticks(x,summary.subject);ax.legend(fontsize=9)
    fig.tight_layout();save(fig,'ica_training_rejection')
    fig,ax=plt.subplots(figsize=(10,3.2));bottom=np.zeros(len(summary))
    for key,label,color in [('eye_blink','Eye','#639c82'),('muscle_artifact','Muscle','#c89156'),('heart_beat','Heart','#a9789e'),('channel_noise','Channel noise','#999999')]:
        values=summary['proposed_'+key].to_numpy();ax.bar(x,values,bottom=bottom,label=label,color=color);bottom+=values
    assert int(bottom.sum())==178
    ax.set(ylabel='Automatic component exclusions',xlabel='Subject',ylim=(0,34));ax.set_xticks(x,summary.subject);ax.legend(ncol=4,fontsize=9)
    fig.tight_layout();save(fig,'ica_component_proposals')
    fig,axes=plt.subplots(1,2,figsize=(10.8,3.3));records=[]
    for ax,(sub,start,ch) in zip(axes,[('005',300,'Fp2'),('014',580,'AF8')]):
        if args.reference:
            sample=np.load(base/f'plot_data/sub-{sub}_blink.npz')
            a=sample['before_uv'];b=sample['after_uv'];times=sample['times']
        else:
            prefix=ROOT/f'artifacts/postica/sub-{sub}/eeg/sub-{sub}_task-ContinuousVideoGamePlay_run-01'
            before=mne.io.read_raw_fif(str(prefix)+'_proc-filt_raw.fif',preload=False,verbose='error').crop(start,start+3.998).load_data()
            after=mne.io.read_raw_fif(str(prefix)+'_proc-clean_raw.fif',preload=False,verbose='error').crop(start,start+3.998).load_data()
            before.set_eeg_reference('average',projection=True,verbose='error').apply_proj(verbose='error')
            assert ch not in before.info['bads']
            a=before.get_data(picks=[ch])[0]*1e6;b=after.get_data(picks=[ch])[0]*1e6;times=before.times+start
        ax.plot(times,a-np.median(a),label='Before ICA',color='#147d92',lw=1)
        ax.plot(times,b-np.median(b),label='After ICA',color='#b75538',lw=1)
        ax.set(title=f'Subject {sub} · {ch}',ylabel='EEG (µV)',xlabel='Time in recording (s)');ax.legend(fontsize=9)
        records.append(dict(subject=sub,start_s=start,channel=ch,before_ptp_uv=float(np.ptp(a)),after_ptp_uv=float(np.ptp(b))))
    fig.tight_layout();save(fig,'ica_blink_attenuation')
    (AUDIT/'figure_provenance.json').write_text(json.dumps(dict(median_rejection_percent=median,automatic_exclusions=178,channel_flags=68,blink_examples=records,
        figure_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for name in ['bad_channel_topomap','ica_training_rejection','ica_component_proposals','ica_blink_attenuation'] for p in [OUT/f'{name}.svg']}),indent=2)+'\n')
    frame=pd.read_csv(base/'postica/results_summary.tsv',sep='\t',dtype={'subject':str})
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True)
    axes[0].bar(frame.subject,frame.rejected_percent,color='#b75538')
    axes[0].set_ylabel('Analysis epochs rejected (%)');axes[0].tick_params(labelbottom=True)
    axes[1].bar(frame.subject,frame.repaired_percent,color='#147d92')
    axes[1].set(ylabel='Retained epochs repaired (%)',xlabel='Subject')
    fig.suptitle('Post-ICA AutoReject on event-related epochs');fig.tight_layout();save(fig,'postica_rejection')
    print('Created five report figures from '+('saved reference inputs.' if args.reference else 'recomputed results.'))

if __name__=='__main__':main()
