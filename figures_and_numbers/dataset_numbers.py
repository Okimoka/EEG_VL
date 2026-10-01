"""Print dataset event counts, durations and recording-boundary exceptions."""
import csv,json,statistics
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]

def rows(path,delimiter='\t'):
    with path.open() as f:return list(csv.DictReader(f,delimiter=delimiter))

def event_code(row):
    value=row['value'].replace('S','').strip()
    return int(value) if value.isdigit() else None

def main():
    dataset=ROOT/'v1.0.0'
    event_tables={};durations={'01':[],'02':[]};counts=[];game=[]
    old_ids=list(range(601,614))+list(range(615,619))
    for number,old in enumerate(old_ids,1):
        sub=f'{number:03d}';folder=dataset/f'sub-{sub}/eeg'
        for run in ('01','02'):
            events=rows(next(folder.glob(f'*run-{run}_events.tsv')));event_tables[sub,run]=events
            channels=rows(next(folder.glob(f'*run-{run}_channels.tsv')))
            fdt=next(folder.glob(f'*run-{run}_eeg.fdt'));sf=json.loads(fdt.with_suffix('.json').read_text())['SamplingFrequency']
            duration=fdt.stat().st_size/(4*len(channels)*sf);durations[run].append(duration)
        c=Counter(r['trial_type'] for run in ('01','02') for r in event_tables[sub,run])
        counts.append(dict(subject=sub,standard=c['ODDBALL STANDARD'],rare=c['ODDBALL RARE'],win=c['GAMBLING WIN'],loss=c['GAMBLING LOSS']))
        logpath=next((dataset/'code/Logs').glob(f'Axon_GAME_Log_{old}_*.csv'))
        with logpath.open() as f:
            table=list(csv.reader(f))[2:]
        log=[(float(r[0]),int(float(r[8]))) for r in table if r and r[0] and float(r[8])!=0]
        first_round=101 if sub=='003' else 102 if sub=='008' else 100
        eeg=event_tables[sub,'02'];index=next(i for i,r in enumerate(eeg) if event_code(r)==first_round)
        j=next(i for i,r in enumerate(log) if r[1]==first_round)
        paired=list(zip(eeg[index:index+50],log[j:j+50]))
        offsets=[b[0]-float(a['onset']) for a,b in paired]
        game.append(dict(subject=sub,last_logged_game_event_min=log[-1][0]/60,
            first_round_marker_s=float(eeg[index]['onset']),log_offset_s=statistics.median(offsets)))
    print('Subject  Standard  Rare  Win  Loss')
    for row in counts:print(f"{row['subject']:>7} {row['standard']:9} {row['rare']:5} {row['win']:4} {row['loss']:5}")
    print('\nMean full recording durations (s):', {r:round(statistics.mean(v),2) for r,v in durations.items()})
    print('Missing gambling:',', '.join(r['subject'] for r in counts if r['win']+r['loss']==0))
    print('Rare-event count range:',min(r['rare'] for r in counts),max(r['rare'] for r in counts))
    print('Subject 001 gameplay through last logged event (min):',round(game[0]['last_logged_game_event_min'],2))
    print('Median first-round marker in the other 15 recordings (s):',round(statistics.median(r['first_round_marker_s'] for r in game if r['subject'] not in ('003','008')),2))
    print('Subject 003 missing initial gameplay (s):',round(game[2]['log_offset_s'],3))
    events=event_tables['008','01'];i=next(i for i,r in enumerate(events) if event_code(r)==100);tail=events[i:]
    logpath=next((dataset/'code/Logs').glob('Axon_GAME_Log_608_*.csv'))
    with logpath.open() as f:table=list(csv.reader(f))[2:]
    log=[(float(r[0]),int(float(r[8]))) for r in table if r and r[0] and float(r[8])!=0]
    offset=statistics.median(float(a['onset'])-b[0] for a,b in zip(tail,log))
    covered=durations['01'][7]-offset;gap=game[7]['log_offset_s']-covered
    print(f'Subject 008: {covered:.3f} s of gameplay in run 01; {gap:.3f} s gap before run 02; run 02 begins {game[7]["log_offset_s"]:.3f} s into the session.')
    print('Subject 008 run-01 tail:',dict(Counter(r['trial_type'] for r in tail)))
    print('No oddball/gambling in run 02:',not any(r['trial_type'].startswith(('ODDBALL','GAMBLING')) for (s,run),events in event_tables.items() if run=='02' for r in events))
    print('The paper reports ~31 min mean gameplay; behavioural logs can continue after playing, so their final timestamp is not a mean active-play duration.')

if __name__=='__main__':main()
