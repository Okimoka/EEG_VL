"""Run every installed PyPREP detector on each filtered recording; save evidence.

No channel is removed/interpolated here. Results require manual review.
"""
import argparse
import contextlib
import hashlib
import inspect
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import time

# Each recording is parallelized separately; do not multiply BLAS threads.
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
import mne
import numpy as np
import pyprep
from pyprep import NoisyChannels
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "pyprep"
SEED = 2026
RANSAC_SETTINGS = dict(n_samples=1000, sample_prop=.25, corr_thresh=.75,
                       frac_bad=.4, corr_window_secs=5., channel_wise=False,
                       max_chunk_size=None)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def software_provenance():
    package = Path(inspect.getfile(pyprep)).parent
    sources = {str(p.relative_to(package)): sha256_file(p)
               for p in sorted(package.rglob("*.py"))}
    return dict(versions={name: version(name) for name in ("pyprep", "mne", "numpy", "scipy")},
                pyprep_source_sha256=hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
                pyprep_source_files=sources)


def atomic_json(path, value):
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(finite(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def archive_previous_results():
    """Keep the complete old proposals before any different recipe replaces them."""
    records = list(OUT.glob("sub-*.json"))
    if not records or all(json.loads(p.read_text()).get("settings", {}).get("random_state") == SEED
                          and json.loads(p.read_text()).get("settings", {}).get("n_samples") == RANSAC_SETTINGS["n_samples"]
                          for p in records):
        return None
    previous = ROOT / "artifacts" / "history"
    previous.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = previous / f"pyprep_before_seed-{SEED}_subsets-{RANSAC_SETTINGS['n_samples']}_{stamp}"
    shutil.copytree(OUT, target)
    atomic_json(target / "archive.json", dict(reason="Preserved before the seed-2026, 1000-subset rerun.",
                created_utc=stamp, recording_count=len(records)))
    return target


def finite(value):
    if isinstance(value, dict):
        return {str(k): finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def run_one(path, force=False):
    subject = path.name.split("sub-")[1].split("_")[0]
    run = path.name.split("run-")[1].split("_")[0]
    key = f"sub-{subject}_run-{run}"
    result_path = OUT / f"{key}.json"
    provenance = software_provenance()
    fingerprint = dict(size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                       input_sha256=sha256_file(path),
                       script_sha256=sha256_file(Path(__file__)),
                       software=provenance,
                       random_state=SEED, ransac=RANSAC_SETTINGS)
    diagnostics_path = OUT / f"{key}_diagnostics.npz"
    if result_path.exists() and not force:
        old = json.loads(result_path.read_text())
        if (old.get("fingerprint") == fingerprint and diagnostics_path.exists()
                and old.get("diagnostics_sha256") == sha256_file(diagnostics_path)):
            return key, "cached"
    start_time = time.monotonic()
    with (OUT / f"{key}.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log), threadpool_limits(limits=1):
        mne.set_log_level("WARNING")
        raw = mne.io.read_raw_fif(path, preload=True, verbose="error")
        if raw.info["bads"]:
            raise ValueError("Detector inputs must precede manual decisions; found bad channels in input.")
        # EEG only; EOG is inspected separately, never spatially interpolated.
        eeg = raw.copy().pick("eeg")
        noisy = NoisyChannels(eeg, do_detrend=True, random_state=SEED,
                              ransac=False, correlation=True, reject_by_annotation=None)
        noisy.find_all_bads(ransac=False, channel_wise=False)
        # find_all_bads does not expose n_samples. Keep every non-RANSAC test,
        # then call the same native RANSAC method with its explicit settings.
        noisy.random_state = np.random.RandomState(SEED)
        excluded = set(noisy.bad_by_nan + noisy.bad_by_flat + noisy.bad_by_manual
                       + noisy.bad_by_correlation + noisy.bad_by_deviation + noisy.bad_by_dropout)
        noisy.find_bad_by_ransac(**RANSAC_SETTINGS)
        bads = noisy.get_bads(as_dict=True, verbose=False)
        chs = list(noisy.ch_names_original)
        duration = raw.n_times / raw.info["sfreq"]
        extra = noisy._extra_info  # Version-pinned diagnostics; keep full arrays.
        arrays = {}
        for criterion, entries in extra.items():
            for name, arr in entries.items():
                arrays[f"{criterion}__{name}"] = np.asarray(arr)
        tested = np.array([ch not in excluded for ch in chs], dtype=bool)
        arrays["bad_by_ransac__tested_channels"] = tested
        arrays["channel_names"] = np.asarray(chs)
        ransac_correlations = arrays["bad_by_ransac__ransac_correlations"]
        if not np.isfinite(ransac_correlations[:, tested]).all():
            raise ValueError("Non-finite correlations in channels actually tested by RANSAC.")
        expected_flags = {ch for ch, fraction in zip(chs, np.mean(ransac_correlations < .75, axis=0))
                          if ch not in excluded and fraction > .4}
        if expected_flags != set(bads["bad_by_ransac"]):
            raise ValueError("Saved RANSAC diagnostics do not reproduce the native flag list.")
        temporary_npz = diagnostics_path.with_name(f".{diagnostics_path.stem}.{os.getpid()}.tmp.npz")
        np.savez_compressed(temporary_npz, **arrays)
        temporary_npz.replace(diagnostics_path)
        windows, metrics = [], {}
        specs = [("bad_by_correlation", "max_correlations", 1., .4, True),
                 ("bad_by_dropout", "dropouts", 1., .0, False),
                 ("bad_by_ransac", "ransac_correlations", 5., .75, True)]
        for criterion, name, seconds, threshold, below in specs:
            arr = np.asarray(extra[criterion][name])
            if criterion == "bad_by_ransac":
                arr = arr.T
            metrics[criterion] = {}
            for idx, ch in enumerate(chs):
                # Native PyPREP fills excluded channels with correlation=1.
                # That is a placeholder, not evidence of a successful test.
                if criterion == "bad_by_ransac" and ch in excluded:
                    metrics[criterion][ch] = None
                    continue
                failed = arr[idx] < threshold if below else arr[idx] > threshold
                metrics[criterion][ch] = float(np.mean(failed))
                # Preserve all failing windows, including currently unflagged channels.
                for w in np.flatnonzero(failed):
                    offset = 0  # PyPREP slices [w*win_size:(w+1)*win_size] for all three tests.
                    windows.append(dict(channel=ch, criterion=criterion,
                                        start=float(w * seconds + offset),
                                        stop=float((w + 1) * seconds + offset),
                                        value=float(arr[idx, w]), threshold=threshold,
                                        detail=("Flag requires >40% failing 5-s windows." if seconds == 5
                                                else "Flag requires >1% failing 1-s windows.")))
        global_specs = [("bad_by_deviation", "robust_channel_deviations", 5),
                        ("bad_by_hf_noise", "hf_noise_zscores", 5),
                        ("bad_by_psd", "psd_zscore", 3)]
        for criterion, name, threshold in global_specs:
            vals = extra[criterion][name]
            metrics[criterion] = dict(zip(chs, map(float, vals)))
            for ch in bads[criterion]:
                windows.append(dict(channel=ch, criterion=criterion, start=0., stop=duration,
                                    value=float(vals[chs.index(ch)]), threshold=threshold,
                                    detail=("Whole-record statistic; no unique trigger interval. "
                                            "PSD can also flag 1/f violations (ratios are diagnostic only); "
                                            "The displayed max ABSOLUTE band z-score includes negative low-power outliers, "
                                            "which do not trigger the band rule. A displayed value >3 alone "
                                            "does not reproduce the flag; band flags require POSITIVE z>3."
                                            if criterion == "bad_by_psd" else
                                            "Whole-record statistic; no unique trigger interval.")))
        for criterion in ("bad_by_nan", "bad_by_flat", "bad_by_SNR"):
            for ch in bads[criterion]:
                windows.append(dict(channel=ch, criterion=criterion, start=0., stop=duration,
                                    value=None, threshold=None,
                                    detail="Whole-record flag; SNR is the intersection of HF noise and low correlation."))
        result = dict(schema_version=2, subject=subject, run=run, recording=key,
                      task="ContinuousVideoGamePlay", sfreq=raw.info["sfreq"],
                      duration=duration, ch_names=raw.ch_names, eeg_ch_names=chs,
                      bads=bads, windows=windows, metrics=metrics,
                      summary={"bad_by_psd": "Metric is max ABSOLUTE band z. Negative low-power outliers can exceed 3 without a PSD flag. Final flags require positive band z>3 or a 1/f violation; ratios are diagnostic only."},
                      all_bads=noisy.get_bads(verbose=False),
                      input=str(path.relative_to(ROOT)), fingerprint=fingerprint,
                      software=provenance, diagnostics_sha256=sha256_file(diagnostics_path),
                      ransac_tested_channels=[ch for ch in chs if ch not in excluded],
                      ransac_excluded_channels=[ch for ch in chs if ch in excluded],
                      ransac_exclusion_note="NaN/flat/manual exclusions and channels already flagged by deviation, correlation or dropout are not tested by RANSAC; native diagnostic correlations of 1 for these channels are placeholders.",
                      settings=dict(random_state=SEED, ransac=True, **RANSAC_SETTINGS,
                                    do_detrend=True, correlation=True, reject_by_annotation=None,
                                    detector_input="0.1–100 Hz, 60 Hz notch, original reference; internal PyPREP 1 Hz detrending",
                                    resampling=None, manual_decisions_applied=False),
                      eog_qc={ch: dict(std_uV=float(np.std(raw.get_data(picks=[ch])) * 1e6),
                                       peak_to_peak_uV=float(np.ptp(raw.get_data(picks=[ch])) * 1e6))
                              for ch in ("VEOG", "HEOG")},
                      seconds_elapsed=round(time.monotonic()-start_time, 2))
        atomic_json(result_path, result)
        return key, f"{len(result['all_bads'])} flags, {result['seconds_elapsed']} s"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--subject", help="Optional subject for debugging, e.g. 001")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted((ROOT / "artifacts" / "mne-bids-pipeline").glob("sub-*/eeg/*_proc-filt_raw.fif"))
    if args.subject:
        files = [p for p in files if p.name.startswith(f"sub-{args.subject}_")]
    expected = 2 if args.subject else 34
    if len(files) != expected:
        raise RuntimeError(f"Expected {expected} filtered recordings, found {len(files)}. Complete pipeline first.")
    if args.jobs < 1:
        p.error("--jobs must be positive")
    archived = archive_previous_results()
    if archived:
        print(f"Archived previous results: {archived.relative_to(ROOT)}", flush=True)
    print(f"PyPREP: seed={SEED}, n_samples={RANSAC_SETTINGS['n_samples']}, {len(files)} recordings, {args.jobs} jobs", flush=True)
    errors = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, f, args.force): f for f in files}
        for future in as_completed(futures):
            try:
                print(*future.result(), flush=True)
            except Exception as exc:
                errors[str(futures[future])] = repr(exc)
                print(f"FAILED {futures[future].name}: {exc}", flush=True)
    atomic_json(OUT / "errors.json", errors)
    if errors:
        raise SystemExit(f"{len(errors)} recordings failed; inspect artifacts/pyprep/errors.json")


if __name__ == "__main__":
    main()
