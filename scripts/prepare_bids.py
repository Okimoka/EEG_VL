"""Create a non-destructive, MNE-readable view of the supplied BIDS dataset.

The default mode repairs channel types and the audited BESA electrode geometry.
--analysis-ready also exports accepted subject-union bad channels and independent
analysis event sidecars (+40 ms once, with close wall/star repeats renamed).
Signal files remain links to the original, unmodified recordings.
"""
import argparse
import csv
import hashlib
import json
import io
import os
from pathlib import Path

from _safe_views import link_source, real_directory, write_independent, remove_linked_locks
from _montage import corrected_coordsystem, corrected_electrodes_tsv, montage_provenance
from _prepared_metadata import (CHANNEL_DESCRIPTION, ANALYSIS_LABELS, EXCLUDED_PREFIX,
                                prepare_events, events_description, prepare_channels,
                                recording_key, verified_channel_policy)

ROOT = Path(__file__).resolve().parents[1]


def _rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _tsv(rows):
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return text.getvalue()


def prepare(source, destination, *, analysis_ready=False,
            policy_path=ROOT / "artifacts/review/channel_policy.json",
            proposal_root=ROOT / "artifacts/pyprep", manifest_path=None):
    source = source.resolve(strict=True)
    destination = Path(os.path.abspath(destination))
    resolved_destination = destination.resolve()
    if (source == resolved_destination or source in resolved_destination.parents
            or resolved_destination in source.parents):
        raise ValueError("Source and prepared dataset directories must not overlap")
    if analysis_ready and resolved_destination in (
            ROOT / "prepared_bids", ROOT / "reviewed_bids"):
        raise ValueError("Use a separate --destination for analysis-ready preparation; preserve prior pipeline inputs")
    if not analysis_ready and destination.exists():
        for path in destination.rglob("*_events.tsv"):
            with path.open(newline="") as handle:
                if "source_event_row" in next(csv.reader(handle, delimiter="\t"), []):
                    raise ValueError("Destination is already analysis-ready; keep --analysis-ready or use a different destination")
    # Validate every modified input before changing a destination sidecar.
    prepared_text = {p: corrected_electrodes_tsv(p)
                     for p in source.rglob("*_electrodes.tsv")}
    for path in source.rglob("*_coordsystem.json"):
        prepared_text[path] = json.dumps(corrected_coordsystem(json.loads(path.read_text())), indent=2) + "\n"
    channel_paths = sorted(source.rglob("*_channels.tsv"))
    event_summary = {}
    policy = None
    if analysis_ready:
        keys = [recording_key(path) for path in channel_paths]
        policy = verified_channel_policy(policy_path, proposal_root, keys)
        for path in sorted(source.rglob("*_events.tsv")):
            eeg_json = path.with_name(path.name.replace("_events.tsv", "_eeg.json"))
            sfreq = json.loads(eeg_json.read_text())["SamplingFrequency"]
            rows, summary = prepare_events(_rows(path), sfreq)
            prepared_text[path] = _tsv(rows)
            metadata_path = path.with_suffix(".json")
            if not metadata_path.exists():
                raise ValueError(f"Expected per-recording event metadata: {metadata_path}")
            prepared_text[metadata_path] = json.dumps(events_description(
                json.loads(metadata_path.read_text())), indent=2) + "\n"
            event_summary[recording_key(path)] = summary
    for path in channel_paths:
        rows = _rows(path)
        if analysis_ready:
            rows = prepare_channels(rows, recording_key(path), policy)
        else:
            for row in rows:
                row.update(type="EOG" if row["name"] in ("VEOG", "HEOG") else "EEG",
                           units="uV")  # EEGLAB storage convention; EOG gain unresolved.
        prepared_text[path] = _tsv(rows)
    for path in source.rglob("*_eeg.json"):
        obj = json.loads(path.read_text())
        obj.update(EEGChannelCount=63, EOGChannelCount=2, MiscChannelCount=0,
                   DigitizedLandmarks=False, DigitizedHeadPoints=False)
        prepared_text[path] = json.dumps(obj, indent=2) + "\n"
    real_directory(destination)
    remove_linked_locks(destination)
    manifest = []
    for original in sorted(source.rglob("*")):
        if not original.is_file() or original.name.endswith(".lock"):
            continue
        target = destination / original.relative_to(source)
        real_directory(target.parent)
        if original in prepared_text:
            write_independent(target, prepared_text[original])
            manifest.append(dict(file=str(original.relative_to(source)),
                                 source_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
                                 prepared_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
        else:
            link_source(original, target)
    generated = []
    if analysis_ready:
        for path in channel_paths:
            target = destination / path.relative_to(source).with_suffix(".json")
            obj = json.loads(path.with_suffix(".json").read_text()) if path.with_suffix(".json").exists() else {}
            obj.update(CHANNEL_DESCRIPTION)
            write_independent(target, json.dumps(obj, indent=2) + "\n")
            generated.append(dict(file=str(target.relative_to(destination)),
                                  prepared_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    if manifest_path is None:
        manifest_path = (ROOT / "artifacts/preparation" / f"{destination.name}.json"
                         if analysis_ready else ROOT / "artifacts/preparation_manifest.json")
    receipt = dict(source=str(source), destination=str(destination),
                   mode="analysis-ready" if analysis_ready else "geometry-and-types",
                   montage_correction=montage_provenance(), modifications=manifest)
    if generated:
        receipt["generated_sidecars"] = generated
    if analysis_ready:
        receipt.update(
            channel_policy_file=str(Path(policy_path).resolve()),
            channel_policy_sha256=hashlib.sha256(Path(policy_path).read_bytes()).hexdigest(),
            proposal_sha256=policy["proposal_sha256"],
            original_per_recording_bads=policy["per_recording_bads"],
            prepared_subject_union_bads=policy["subject_union_bads"],
            event_rules=dict(shift_seconds=0.040, shifted_trial_types=sorted(ANALYSIS_LABELS),
                             unchanged_timing="All other labels, including STATUS/boundary and task-start markers",
                             repeated_trial_types=["PLAYER_CRASH_WALL", "COLLECT_STAR"],
                             minimum_same_type_separation_seconds=0.500,
                             comparison="Original samples; compare with last retained event of the same type",
                             excluded_prefix=EXCLUDED_PREFIX,
                             value_column="Original trigger values unchanged; analysis selection uses prepared trial_type",
                             source_rows="Preserved with one-based source_event_row; no EEG sample changes",
                             ica_effect="None for native task_is_rest fixed-window selection"),
            events=event_summary)
    write_independent(manifest_path, json.dumps(receipt, indent=2) + "\n")
    print(f"Prepared {destination}; wrote {len(manifest)} independent metadata files; "
          + ("analysis events and shared bad-channel sets ready." if analysis_ready else "events unchanged."))
    print(f"Provenance: {manifest_path}")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "v1.0.0")
    parser.add_argument("--destination", type=Path, default=ROOT / "prepared_bids")
    parser.add_argument("--analysis-ready", action="store_true",
                        help="Export analysis events and accepted common bad channels; requires a separate destination")
    args = parser.parse_args()
    prepare(args.source, args.destination, analysis_ready=args.analysis_ready)
