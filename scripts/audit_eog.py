#!/usr/bin/env python3
"""Small diagnostic audit; does not modify any recording or apply a gain correction.

Units reproduce MNE's EEGLAB assumption (native values * 1e-6), not a verified
calibration for the bipolar auxiliary channels. Raw quantiles use every 50th
sample; filtering uses a separate middle 120-second segment at native 500 Hz.
"""
from pathlib import Path
import csv,json,numpy as np
from scipy.signal import butter,sosfiltfilt
root=Path(__file__).resolve().parents[1]; rows=[]
for p in sorted((root/'v1.0.0').glob('sub-*/eeg/*_eeg.fdt')):
 channel_file=p.with_name(p.name.replace('_eeg.fdt','_channels.tsv'))
 channels=[r['name'] for r in csv.DictReader(channel_file.open(),delimiter='\t')]
 sf=float(json.loads(p.with_suffix('.json').read_text())['SamplingFrequency'])
 a=np.memmap(p,mode='r',dtype='<f4').reshape(-1,len(channels))
 indices=[channels.index(n) for n in ['Fp1','Pz','VEOG','HEOG']]
 sampled=np.asarray(a[::50,indices],dtype=float)*1e-6
 nwin=min(len(a),int(120*sf));start=(len(a)-nwin)//2
 segment=np.asarray(a[start:start+nwin,indices],dtype=float).T*1e-6
 fil=sosfiltfilt(butter(4,[.5,15],btype='bandpass',fs=sf,output='sos'),segment)
 # Discard 5-s window edges for finite-segment filter diagnostic.
 fil=fil[:,int(5*sf):-int(5*sf)]
 fp1_minus_pz=fil[0]-fil[1]
 for i,ch in enumerate(['Fp1','Pz','VEOG','HEOG']):
  x=sampled[:,i];q=np.quantile(x,[.01,.5,.99]);qf=np.quantile(fil[i],[.01,.99])
  rows.append(dict(recording=p.stem.replace('_eeg',''),channel=ch,n_samples=len(a),sampled_every_n=50,unit='V_under_MNE_EEGLAB_assumption',median=q[1],p01=q[0],p99=q[2],center_segment_start_s=start/sf,center_segment_length_s=nwin/sf,bandpass_0p5_15_p99_minus_p01=qf[1]-qf[0],bandpass_corr_fp1_minus_pz=np.corrcoef(fil[i],fp1_minus_pz)[0,1],all_sampled_finite=bool(np.isfinite(x).all())))
 print(p.stem,'V medians',np.round(np.median(sampled,axis=0),4),'bandpassed EOG p99-p01',np.round([rows[-2]['bandpass_0p5_15_p99_minus_p01'],rows[-1]['bandpass_0p5_15_p99_minus_p01']],4),flush=True)
output=root/'artifacts/qc'
output.mkdir(parents=True,exist_ok=True)
with (output/'eog_channel_audit.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
