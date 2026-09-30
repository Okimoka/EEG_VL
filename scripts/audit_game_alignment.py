from pathlib import Path
import pandas as pd,numpy as np,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'artifacts/report_notes_revision'
OUT.mkdir(parents=True,exist_ok=True)
results=[]
for subject,old,start in [('003',603,101),('008',608,102)]:
 log=pd.read_csv(next((ROOT/'v1.0.0/code/Logs').glob(f'Axon_GAME_Log_{old}_*.csv')),skiprows=[1],usecols=[0,8]);log.columns=['seconds','code'];log=log[log.code.ne(0)].reset_index(drop=True)
 def load(run):
  e=pd.read_csv(next((ROOT/f'v1.0.0/sub-{subject}/eeg').glob(f'*run-{run}_events.tsv')),sep='\t')
  e['code']=pd.to_numeric(e.value.astype(str).str.replace('S','',regex=False).str.strip(),errors='coerce');return e
 eeg=load('02');eeg=eeg.iloc[np.flatnonzero(eeg.code.eq(start))[0]:].reset_index(drop=True)
 matching=log.iloc[np.flatnonzero(log.code.eq(start))[0]:].reset_index(drop=True)
 assert np.array_equal(eeg.code.head(50),matching.code.head(50))
 differences=matching.seconds.head(50).to_numpy()-eeg.onset.head(50).to_numpy()
 result=dict(subject=subject,matched_first_50_codes=True,offset_s=float(np.median(differences)),offset_range_s=[float(differences.min()),float(differences.max())])
 if subject=='008':
  first=load('01');first=first.iloc[np.flatnonzero(first.code.eq(100))[0]:].reset_index(drop=True)
  assert np.array_equal(first.code,log.code.head(len(first)))
  diff=first.onset.to_numpy()-log.seconds.head(len(first)).to_numpy()
  recordings=pd.read_csv(ROOT/'artifacts/qc/recordings.tsv',sep='\t',dtype={'subject':str,'run':str})
  duration=float(recordings.loc[recordings.subject.eq(subject)&recordings.run.eq('01'),'duration_s'].iloc[0])
  covered=duration-float(np.median(diff))
  result.update(run1_tail_codes_matched=int(len(first)),run1_tail_recorded_s=covered,gap_between_recordings_s=float(np.median(differences))-covered)
 results.append(result)
(OUT/'game_alignment.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
