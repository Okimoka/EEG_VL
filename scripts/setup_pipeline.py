"""Apply/verify the recorded upstream compatibility change in the pinned environment.

Source: https://github.com/mne-tools/mne-bids-pipeline/pull/1164
No scientific data or configuration is changed.
"""
import argparse

import ast

import copy

from datetime import datetime, timezone

import difflib

import hashlib

from importlib.metadata import distribution

import json

from pathlib import Path

import tempfile

from types import SimpleNamespace

RELATIVE = "mne_bids_pipeline/steps/preprocessing/_06a1_fit_ica.py"

BEFORE_SHA256 = "f5e522fd26c75b4a291b21cd77c33925a987be34526b53ff4aea0f7e33945553"

AFTER_SHA256 = "33af95feae4f3caf548fc231279059ce14bba4314f5d85f77497ae41763387a0"

OLD = '''            if cfg.ica_l_freq is not None or h_freq is not None:
                logger.info(**gen_log_kwargs(message=msg))
                raw.filter(l_freq=cfg.ica_l_freq, h_freq=h_freq, n_jobs=1)
            del nyq, h_freq
'''

NEW = '''            if msg:
                logger.info(**gen_log_kwargs(message=msg))
            del nyq

        if cfg.ica_l_freq is not None or h_freq is not None:
            raw.filter(l_freq=cfg.ica_l_freq, h_freq=h_freq, n_jobs=1)
'''

def digest(data):
    return hashlib.sha256(data).hexdigest()

def exercise_filter_loop(source):
    """Execute the installed loop up to epoch creation against instrumented raws.

    The AST comes from the real package source, including its first-run branch,
    Nyquist handling, and filter call. Only I/O/logging objects are mocked.
    """
    module = ast.parse(source)
    function = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == "run_ica")
    loop = next(n for n in function.body if isinstance(n, ast.For)
                and ast.unparse(n.iter) == "enumerate(zip(cfg.runs, raw_fnames))")
    loop = copy.deepcopy(loop)
    stop = next(i for i, n in enumerate(loop.body) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "event_id" for t in n.targets))
    loop.body = loop.body[:stop]
    code = compile(ast.fix_missing_locations(ast.Module(body=[loop], type_ignores=[])),
                   "<installed ICA prefilter loop>", "exec")
    scenarios = [
        ("two_runs_1_to_100_Hz", 500., 1., 100., True, 500., (1., 100.)),
        ("ICLabel_Nyquist_fallback_both_runs", 128., 1., 100., True, None, (1., None)),
        ("lowpass_only_both_runs", 500., None, 40., False, None, (None, 40.)),
        ("disabled_filter_stays_disabled", 500., None, None, False, None, None),
    ]
    checks = []
    for name, sfreq, low, high, iclabel, resample, expected_filter in scenarios:
        calls, reads, warnings = [], [], []

        class InstrumentedRaw:
            def __init__(self, name):
                self.name = name
                self.info = {"sfreq": sfreq}

            def filter(self, **kwargs):
                calls.append((self.name, kwargs))
                return self

        def read_raw(fname, preload):
            assert preload is True
            reads.append(fname.basename)
            return InstrumentedRaw(fname.basename)

        environment = dict(
            cfg=SimpleNamespace(runs=["01", "02"], ica_l_freq=low, ica_h_freq=high,
                                ica_use_icalabel=iclabel, raw_resample_sfreq=resample),
            raw_fnames=[SimpleNamespace(basename=run) for run in ("01", "02")],
            mne=SimpleNamespace(io=SimpleNamespace(read_raw_fif=read_raw)),
            np=SimpleNamespace(allclose=lambda a, b: a == b),
            logger=SimpleNamespace(info=lambda **kw: None, warning=lambda **kw: warnings.append(kw)),
            gen_log_kwargs=lambda **kw: kw,
        )
        exec(code, environment)
        expected = [] if expected_filter is None else [
            (run, dict(l_freq=expected_filter[0], h_freq=expected_filter[1], n_jobs=1))
            for run in ("01", "02")]
        if reads != ["01", "02"] or calls != expected:
            raise AssertionError(f"{name}: expected {expected}, got {calls}; reads={reads}")
        if name.startswith("ICLabel_") and len(warnings) != 1:
            raise AssertionError("Expected one first-run Nyquist warning")
        checks.append(dict(scenario=name, runs_read=reads, filter_calls=len(calls), passed=True))
    return checks

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    target = Path(distribution('mne-bids-pipeline').locate_file(RELATIVE))
    data = target.read_bytes()
    if digest(data) == BEFORE_SHA256:
        if args.check:
            raise SystemExit('Run this script without --check to finish environment setup.')
        source = data.decode()
        assert source.count(OLD) == 1
        data = source.replace(OLD, NEW).encode()
        assert digest(data) == AFTER_SHA256
        exercise_filter_loop(data.decode())
        target.write_bytes(data)
    if digest(data) != AFTER_SHA256:
        raise SystemExit('Unexpected pipeline source; install the pinned requirements first.')
    checks = exercise_filter_loop(data.decode())
    print(f'Pipeline environment verified; {len(checks)} two-run checks passed.')

if __name__ == '__main__':
    main()
