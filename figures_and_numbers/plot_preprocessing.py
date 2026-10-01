"""Draw the bad-channel/ICA summaries and print their numbers."""
import os
from collections import Counter
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLBACKEND','Agg')
from pathlib import Path
import json
import matplotlib.pyplot as plt
import mne,numpy as np,pandas as pd
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'figures'
INPUT = ROOT/'logs/plot_data'
SPHERE = (0, 0.004320136895757215, 0.04231859617345634, .085)

def save(fig,name):
    fig.savefig(OUT/f'{name}.svg',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

def main():
    mne.set_log_level('ERROR')
    proposals=[json.loads(p.read_text()) for p in sorted((ROOT/'logs/pyprep').glob('sub-*_run-*.json'))]
    total=Counter(ch for result in proposals for ch in result['all_bads'])
    run_counts={run:Counter(ch for result in proposals if result['run']==run for ch in result['all_bads']) for run in ('01','02')}
    electrodes=pd.read_csv(ROOT/'00_prepare_dataset/electrodes.tsv',sep='\t').iloc[:63]
    info=mne.create_info(electrodes.name.tolist(),500,'eeg')
    info.set_montage(mne.channels.make_dig_montage(
        ch_pos=dict(zip(electrodes.name,electrodes[['x','y','z']].to_numpy())),coord_frame='head'))
    counts=pd.DataFrame([dict(channel=ch,all_recordings=total[ch],run_01=run_counts['01'][ch],run_02=run_counts['02'][ch]) for ch in info.ch_names])
    per_recording=[len(result['all_bads']) for result in proposals]
    ventral=sum(total[ch] for ch in ('FT9','TP9','TP10','FT10'))
    stats=dict(flagged_channel_recordings=sum(per_recording),median_flags=float(np.median(per_recording)),
        min_flags=min(per_recording),max_flags=max(per_recording),recordings_without_flags=per_recording.count(0),
        most_frequent=total.most_common(1),ventral_flags=ventral,ventral_cases=34*4,
        other_flags=sum(per_recording)-ventral,other_cases=34*59)
    counts=counts.set_index('channel').loc[info.ch_names]
    vmax=int(counts.all_recordings.max());fig,axes=plt.subplots(1,3,figsize=(11.5,4.0))
    for ax,column,title in zip(axes,['all_recordings','run_01','run_02'],['All recordings (34)','Run 01: oddball / gambling (17)','Run 02: gameplay (17)']):
        names=[f'{name}\n{int(counts.loc[name,column])}' if counts.loc[name,'all_recordings']>=5 else '' for name in info.ch_names]
        im,_=mne.viz.plot_topomap(counts[column].to_numpy(float),info,axes=ax,names=names,
             contours=0,image_interp='linear',cmap='YlOrRd',vlim=(0,vmax),show=False,sphere=SPHERE,extrapolate='head')
        for text in ax.texts:
            if text.get_text():text.set_bbox(dict(facecolor='white',alpha=.8,edgecolor='none',pad=.8))
        ax.set_title(title,fontsize=11)
    fig.subplots_adjust(left=.02,right=.89,top=.87,bottom=.06,wspace=.2)
    cax=fig.add_axes([.92,.2,.018,.55]);bar=fig.colorbar(im,cax=cax);bar.set_label('Recordings with a flag');bar.set_ticks(range(vmax+1))
    save(fig,'bad_channel_topomap')
    summary=pd.read_csv(ROOT/'logs/ica_final_full/results_summary.tsv',sep='\t',dtype={'subject':str})
    x=np.arange(len(summary));median=float(summary.rejected_percent.median())
    fig,ax=plt.subplots(figsize=(10,3.2));ax.bar(x,summary.rejected_percent,color='#597c9c')
    ax.axhline(median,color='black',linestyle=':',linewidth=1,label=f'Median over all subjects ({median:.1f}%)')
    ax.set(ylabel='ICA fitting windows rejected (%)',xlabel='Subject',ylim=(0,55));ax.set_xticks(x,summary.subject);ax.legend(fontsize=9)
    fig.tight_layout();save(fig,'ica_training_rejection')
    fig,ax=plt.subplots(figsize=(10,3.2));bottom=np.zeros(len(summary))
    for key,label,color in [('eye_blink','Eye','#639c82'),('muscle_artifact','Muscle','#c89156'),('heart_beat','Heart','#a9789e'),('channel_noise','Channel noise','#999999')]:
        values=summary['proposed_'+key].to_numpy();ax.bar(x,values,bottom=bottom,label=label,color=color);bottom+=values
    ax.set(ylabel='Automatic component exclusions',xlabel='Subject',ylim=(0,34));ax.set_xticks(x,summary.subject);ax.legend(ncol=4,fontsize=9)
    fig.tight_layout();save(fig,'ica_component_proposals')
    fig,axes=plt.subplots(1,2,figsize=(10.8,3.3));records=[]
    for ax,(sub,start,ch) in zip(axes,[('005',300,'Fp2'),('014',580,'AF8')]):
        extract=np.load(INPUT/f'sub-{sub}_blink.npz')
        a=extract['before_uv'];b=extract['after_uv'];times=extract['times']
        ax.plot(times,a-np.median(a),label='Before ICA',color='#147d92',lw=1)
        ax.plot(times,b-np.median(b),label='After ICA',color='#b75538',lw=1)
        ax.set(title=f'Subject {sub} · {ch}',ylabel='EEG (µV)',xlabel='Time in recording (s)');ax.legend(fontsize=9)
        records.append(dict(subject=sub,start_s=start,channel=ch,before_ptp_uv=float(np.ptp(a)),after_ptp_uv=float(np.ptp(b))))
    fig.tight_layout();save(fig,'ica_blink_attenuation')
    stats.update(ica_candidate_windows=int(summary.candidate_windows.sum()),ica_retained_windows=int(summary.retained_windows.sum()),
        ica_components=int(summary.components.sum()),automatic_exclusions=int(summary.proposed_exclusions.sum()),
        reviewed_exclusions=int(summary.current_status_bad.sum()),median_reviewed_exclusions=float(summary.current_status_bad.median()),
        blink_examples=records)
    print(json.dumps(stats,indent=2))
    (ROOT/'logs/preprocessing_numbers.json').write_text(json.dumps(stats,indent=2)+'\n')


if __name__=='__main__':main()
