"""Invoke the ordinary MNE-BIDS-Pipeline CLI with thread limits and a saved log.

Usage: python scripts/run_pipeline.py --log artifacts/NAME/native_pipeline.log
       --config CONFIG.py --steps STEPS [other native CLI options]
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__,add_help=False)
    parser.add_argument('--log',required=True,type=Path)
    args,native=parser.parse_known_args()
    if not native:parser.error('Supply --config and --steps for the native pipeline.')
    environment=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        environment[name]='1'
    environment.update(MPLBACKEND='Agg',PYTHONUNBUFFERED='1')
    path=ROOT/args.log;path.parent.mkdir(parents=True,exist_ok=True)
    command=[sys.executable,'-c','from mne_bids_pipeline._main import main; main()',*native]
    with path.open('a') as log:
        process=subprocess.Popen(command,cwd=ROOT,env=environment,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        try:
            for line in process.stdout:
                print(line,end='',flush=True);log.write(line);log.flush()
            status=process.wait()
        except BaseException:
            process.terminate();process.wait();raise
        message=f'\npipeline exit status: {status}\n';log.write(message);print(message,end='')
    raise SystemExit(status)

if __name__=='__main__':main()
