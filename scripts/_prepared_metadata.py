"""Explicit analysis metadata for the independent, pipeline-ready BIDS view."""
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import json
from pathlib import Path
import re

ANALYSIS_LABELS = frozenset((
    "ODDBALL STANDARD", "ODDBALL RARE", "GAMBLING WIN", "GAMBLING LOSS",
    "PLAYER_CRASH_WALL", "COLLECT_STAR", "MISSILE_HIT_ENEMY",
    "PLAYER_CRASH_ENEMY", "COLLECT_AMMO", "SHOOT_BUTTON",
))
REPEATED_LABELS = frozenset(("PLAYER_CRASH_WALL", "COLLECT_STAR"))
EXCLUDED_PREFIX = "EXCLUDED_REPEAT__"
EVENT_COLUMNS = ("source_event_row", "source_onset", "source_sample", "source_trial_type",
                 "onset_shift_seconds", "keep_same_type_500ms", "exclusion_reason")


def recording_key(path):
    match = re.match(r"(sub-[^_]+)_.*_run-([^_]+)_", Path(path).name)
    if match is None:
        raise ValueError(f"Missing subject/run in {path}")
    return f"{match[1]}_run-{match[2]}"


def prepare_events(rows, sfreq):
    """Keep every source row; shift analysis anchors once and relabel repeats.

    The input must be the original event table. Timing comparisons use original
    sample units; each event type compares against its last retained event.
    Acquisition/status/start markers retain their original timing.
    """
    frequency = Decimal(str(sfreq))
    if frequency <= 0 or frequency % 2:
        raise ValueError("The 500 ms threshold must be an integer number of samples")
    last, result = {}, []
    previous = None
    for number, row in enumerate(rows, 1):
        if any(field in row for field in EVENT_COLUMNS):
            raise ValueError("Events are already prepared; use the original source dataset")
        label = row["trial_type"]
        if label.startswith(EXCLUDED_PREFIX):
            raise ValueError("Input already contains excluded-repeat labels")
        onset = Decimal(str(row["onset"]))
        if not onset.is_finite():
            raise ValueError("Event onsets must be finite")
        if previous is not None and onset < previous:
            raise ValueError("Source event rows must be chronological")
        previous = onset
        sample = int((onset * frequency).to_integral_value(rounding=ROUND_HALF_EVEN))
        keep = True
        if label in REPEATED_LABELS:
            keep = label not in last or sample - last[label] >= int(frequency / 2)
            if keep:
                last[label] = sample
        shift = Decimal("0.040") if label in ANALYSIS_LABELS else Decimal(0)
        out = dict(row)
        out.update(source_event_row=str(number), source_onset=str(row["onset"]),
                   source_sample=str(row.get("sample", "n/a")), source_trial_type=label,
                   onset_shift_seconds=f"{shift:.3f}", keep_same_type_500ms=str(keep).lower(),
                   exclusion_reason="n/a" if keep else "same-type event <500 ms after last retained")
        out["onset"] = f"{onset + shift:.10f}"
        if row.get("sample", "n/a") != "n/a":
            source_sample = Decimal(str(row["sample"]))
            if source_sample != source_sample.to_integral_value():
                raise ValueError("Source sample indices must be integers")
            out["sample"] = str(int(source_sample + shift * frequency))
        if not keep:
            out["trial_type"] = EXCLUDED_PREFIX + label
        result.append(out)
    counts = Counter(row["source_trial_type"] for row in result)
    excluded = Counter(row["source_trial_type"] for row in result
                       if row["keep_same_type_500ms"] == "false")
    return result, {"source_rows": len(result), "shifted_analysis_rows": sum(
        row["source_trial_type"] in ANALYSIS_LABELS for row in result),
        "source_counts": dict(sorted(counts.items())), "excluded_repeated_counts": dict(excluded)}


def events_description(original):
    description = dict(original)
    description["onset"] = dict(Description="Prepared onset: original onset plus 0.040 s for the ten analysis event types; acquisition/status/start markers unchanged", Units="s")
    description["trial_type"] = dict(Description="Original trial type for eligible events; EXCLUDED_REPEAT__ prefix marks wall/star anchors excluded by the retained-first same-type 500 ms rule. Match analysis conditions by their original exact labels.")
    description.update({
        "source_event_row": dict(Description="One-based data row in the original events.tsv (header not counted)"),
        "source_onset": dict(Description="Unmodified original onset", Units="s"),
        "source_sample": dict(Description="Unmodified original sample field; n/a if not supplied"),
        "source_trial_type": dict(Description="Unmodified original trial_type"),
        "onset_shift_seconds": dict(Description="Correction already applied to onset and any numeric sample index; do not apply again", Units="s"),
        "keep_same_type_500ms": dict(Description="False only for wall/star events less than 500 ms after the last retained same-type event, measured on the original sampling grid", Levels={"true": "Retain this anchor", "false": "Excluded repeated anchor; EEG and event row preserved"}),
        "exclusion_reason": dict(Description="Reason for renaming the analysis event, or n/a"),
    })
    return description


def verified_channel_policy(policy_path, proposal_root, recording_keys):
    """Bind accepted decisions to the actual detector files before exporting."""
    policy_path = Path(policy_path)
    policy = json.loads(policy_path.read_text())
    if policy["random_state"] != 2026 or policy["n_samples"] != 1000:
        raise ValueError("Expected the accepted seed-2026, 1000-subset PyPREP policy")
    if set(recording_keys) != set(policy["per_recording_bads"]):
        raise ValueError("Source recordings differ from the accepted channel-policy recordings")
    unions = defaultdict(set)
    for key in recording_keys:
        path = Path(proposal_root) / f"{key}.json"
        if hashlib.sha256(path.read_bytes()).hexdigest() != policy["proposal_sha256"][key]:
            raise ValueError(f"{key}: detector hash differs from the accepted policy")
        proposal = json.loads(path.read_text())
        if proposal["settings"]["random_state"] != 2026 or proposal["settings"]["n_samples"] != 1000:
            raise ValueError(f"{key}: unexpected detector settings")
        if set(proposal["all_bads"]) != set(policy["per_recording_bads"][key]):
            raise ValueError(f"{key}: accepted decisions differ from detector flags")
        unions[key.split("_")[0].removeprefix("sub-")].update(proposal["all_bads"])
    if {key: sorted(value) for key, value in unions.items()} != policy["subject_union_bads"]:
        raise ValueError("Accepted subject union is inconsistent with per-run flags")
    return policy


def prepare_channels(rows, key, policy):
    subject = key.split("_")[0].removeprefix("sub-")
    run_bads = set(policy["per_recording_bads"][key])
    union_bads = set(policy["subject_union_bads"][subject])
    names = {row["name"] for row in rows}
    if not union_bads <= names or union_bads & {"VEOG", "HEOG"}:
        raise ValueError(f"{key}: invalid accepted EEG bad-channel set")
    result = []
    for row in rows:
        row = dict(row)
        name = row["name"]
        row.update(type="EOG" if name in ("VEOG", "HEOG") else "EEG", units="uV",
                   status="bad" if name in union_bads else "good",
                   status_description=("Accepted PyPREP subject union across both runs (seed 2026; RANSAC n_samples=1000)" if name in union_bads else "n/a"),
                   pyprep_bad_in_this_run=str(name in run_bads).lower())
        result.append(row)
    return result


CHANNEL_DESCRIPTION = {
    "pyprep_bad_in_this_run": {"Description": "Accepted original per-run PyPREP flag before taking the subject union; the standard status column uses the same union in both runs", "Levels": {"true": "Flagged in this recording", "false": "Not flagged in this recording"}},
}
