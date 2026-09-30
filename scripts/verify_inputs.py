"""Verify shipped reference files and the exact raw dataset used in this project."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def verify(root,manifest):
    for name,expected in manifest.items():
        path=root/name
        if not path.is_file():raise SystemExit(f'Missing input: {path}')
        with path.open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
        if actual!=expected:raise SystemExit(f'Input differs from the recorded project: {path}')
    print(f'Verified {len(manifest)} files in {root}')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,help='Original dataset root (hashes all signal and sidecar files)')
    args=p.parse_args()
    verify(ROOT,json.loads((ROOT/'reference/file_manifest.json').read_text()))
    if args.data:verify(args.data,json.loads((ROOT/'reference/source_manifest.json').read_text()))

if __name__=='__main__':main()
