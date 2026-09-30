"""Export the adopted all-PyPREP-flags policy; optionally restore saved proposals."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', action='store_true', help='Use the supplied accepted detector JSONs, without rerunning detection')
    args = parser.parse_args()
    proposals = ROOT/'artifacts/pyprep'
    proposals.mkdir(parents=True, exist_ok=True)
    if args.reference:
        for source in (ROOT/'reference/pyprep').glob('sub-*_run-*.json'):
            target = proposals/source.name
            if target.exists() and target.read_bytes() != source.read_bytes():
                raise SystemExit(f'Refusing to overwrite new detector results: {target}')
            shutil.copy2(source, target)
    files = sorted(proposals.glob('sub-*_run-*.json'))
    expected = {f'sub-{s:03d}_run-{r:02d}' for s in range(1,18) for r in (1,2)}
    if {p.stem for p in files} != expected:
        raise SystemExit('Need all 34 PyPREP results.')
    decisions, hashes, unions = {}, {}, defaultdict(set)
    counts = {run: Counter() for run in ('01','02')}
    for path in files:
        result = json.loads(path.read_text())
        assert result['settings']['random_state'] == 2026 and result['settings']['n_samples'] == 1000
        decisions[path.stem] = result['all_bads']
        hashes[path.stem] = hashlib.sha256(path.read_bytes()).hexdigest()
        subject, run = path.stem.removeprefix('sub-').split('_run-')
        unions[subject].update(result['all_bads'])
        counts[run].update(result['all_bads'])
    policy = dict(policy='Accept all PyPREP threshold-based decisions; no manual overrides',
        random_state=2026, n_samples=1000, proposal_sha256=hashes,
        per_recording_bads=decisions, subject_union_bads={s: sorted(ch) for s,ch in unions.items()},
        total_flagged_channel_recordings=sum(map(len,decisions.values())))
    target = ROOT/'artifacts/review/channel_policy.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(policy,indent=2)+'\n')
    qc=ROOT/'artifacts/qc';qc.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([dict(channel=ch,all_recordings=counts['01'][ch]+counts['02'][ch],run_01=counts['01'][ch],run_02=counts['02'][ch])
        for ch in result['eeg_ch_names']]).to_csv(qc/'bad_channel_frequency.tsv',sep='\t',index=False)
    recorded=json.loads((ROOT/'reference/channel_policy.json').read_text())
    differences={key:dict(recorded=recorded['per_recording_bads'][key],new=bad)
        for key,bad in decisions.items() if set(bad)!=set(recorded['per_recording_bads'][key])}
    (ROOT/'artifacts/review/differences_from_reference.json').write_text(json.dumps(differences,indent=2)+'\n')
    print(f"Accepted {policy['total_flagged_channel_recordings']} flags; {len(differences)} recordings differ from reference.")
    if differences:
        print('Saved reference ICA models cannot be assumed compatible; inspect the differences before continuing.')

if __name__=='__main__':main()
