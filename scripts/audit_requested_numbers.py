from pathlib import Path
import pandas as pd,numpy as np,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'artifacts/report_notes_revision'
OUT.mkdir(parents=True,exist_ok=True)
old_ids=[601,602,603,604,605,606,607,608,609,610,611,612,613,615,616,617,618]
rows=[];eyes=[];durations=[]
for number,old in enumerate(old_ids,1):
 sub=f'{number:03d}'
 events=[]
 for run in ['01','02']:
  p=next((ROOT/f'v1.0.0/sub-{sub}/eeg').glob(f'*run-{run}_events.tsv'))
  e=pd.read_csv(p,sep='\t');e['run']=run;events.append(e)
  if run=='01':
   for task in ['ODDBALL','GAMBLING']:
    tt=e[e.trial_type.str.startswith(task)]
    if len(tt):durations.append(dict(subject=sub,task=task,count=len(tt),first_s=float(tt.onset.min()),last_s=float(tt.onset.max()),span_min=float((tt.onset.max()-tt.onset.min())/60)))
 game_path=next((ROOT/'v1.0.0/code/Logs').glob(f'Axon_GAME_Log_{old}_*.csv'))
 log=pd.read_csv(game_path,skiprows=[1],usecols=[0,8]);log.columns=['seconds','code']
 starts=log[log.code.eq(100)]
 e=events[1];codes=pd.to_numeric(e.value.astype(str).str.replace('S','',regex=False).str.strip(),errors='coerce')
 first_round=101 if old==603 else 102 if old==608 else 100
 eeg_start=float(e.loc[codes.eq(first_round),'onset'].iloc[0])
 log_round=float(log.loc[log.code.eq(first_round),'seconds'].iloc[0])
 row=dict(subject=sub,original_id=old,log_duration_min=float((log.seconds.max()-log.seconds.min())/60),game_event_span_min=float((log.loc[log.code.ne(0),'seconds'].max()-log.loc[log.code.ne(0),'seconds'].min())/60),
  first_round_in_run2=first_round,eeg_round_onset_s=eeg_start,log_round_onset_s=log_round,
  implied_log_time_at_eeg_start_s=log_round-eeg_start)
 rows.append(row)
 table=pd.read_csv(next((ROOT/f'artifacts/ica_native_cohort/sub-{sub}/eeg').glob('*components.tsv')),sep='\t')
 eye=table.status_description.str.contains('Auto-detected eye blink',na=False)
 eyes.append(dict(subject=sub,automatic_eye=int(eye.sum()),first10_eye=int((eye & (table.component<10)).sum()),total_excluded=int(table.status.eq('bad').sum())))
pd.DataFrame(rows).to_csv(OUT/'game_timing.tsv',sep='\t',index=False)
pd.DataFrame(durations).to_csv(OUT/'exemplar_event_spans.tsv',sep='\t',index=False)
pd.DataFrame(eyes).to_csv(OUT/'component_counts.tsv',sep='\t',index=False)
r=pd.DataFrame(rows);ee=pd.DataFrame(eyes);d=pd.DataFrame(durations)
f=pd.read_csv(ROOT/'artifacts/ica_native_cohort/results_summary.tsv',sep='\t')
summary=dict(game_log_duration_mean_min=float(r.log_duration_min.mean()),longest_game=r.loc[r.game_event_span_min.idxmax()].to_dict(),longest_log=r.loc[r.log_duration_min.idxmax()].to_dict(),
 median_other_game_start_s=float(r.loc[~r.subject.isin(['003','008']),'eeg_round_onset_s'].median()),
 game_exceptions=r[r.subject.isin(['003','008'])].to_dict('records'),
 mean_event_spans_min=d.groupby('task').span_min.mean().to_dict(),
 total_automatic_eye=int(ee.automatic_eye.sum()),first10_eye_total=int(ee.first10_eye.sum()),first10_eye_mean=float(ee.first10_eye.mean()),
 automatic_eye_mean=float(ee.automatic_eye.mean()),excluded_total=int(ee.total_excluded.sum()),excluded_median=float(ee.total_excluded.median()),
 median_rejection_percent=float(f.rejected_percent.median()),paper_ocular_mean=(13*3+4*2)/17)
(OUT/'requested_numbers.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
