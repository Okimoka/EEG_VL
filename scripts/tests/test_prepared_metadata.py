"""Scientific and source-preservation checks for analysis-ready preparation."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _prepared_metadata import (prepare_events, prepare_channels, verified_channel_policy,
                                EXCLUDED_PREFIX)
from prepare_bids import prepare


def event(onset, label="PLAYER_CRASH_WALL", sample="n/a"):
    return dict(onset=str(onset), duration="n/a", sample=sample,
                trial_type=label, value="S 10")


class PreparedMetadataTests(unittest.TestCase):
    def test_retained_first_chain_exact_boundary_and_independent_types(self):
        rows = [event(0), event(.3), event(.4, "COLLECT_STAR"), event(.5),
                event(.6), event(.9, "COLLECT_STAR"), event(1)]
        prepared, _ = prepare_events(rows, 500)
        self.assertEqual([r["keep_same_type_500ms"] for r in prepared],
                         ["true", "false", "true", "true", "false", "true", "true"])
        self.assertEqual(prepared[1]["trial_type"], EXCLUDED_PREFIX + "PLAYER_CRASH_WALL")
        self.assertEqual(prepared[-1]["onset"], "1.0400000000")
        self.assertEqual(rows[1]["trial_type"], "PLAYER_CRASH_WALL")
        chain, _ = prepare_events([event(0), event(.3), event(.6)], 500)
        self.assertEqual([r["keep_same_type_500ms"] for r in chain], ["true", "false", "true"])

    def test_timing_source_fields_and_no_shift_of_boundary(self):
        rows = [event(0, "STATUS", "0"), event(1, "ODDBALL STANDARD", "500")]
        prepared, _ = prepare_events(rows, 500)
        self.assertEqual(prepared[0]["onset"], "0.0000000000")
        self.assertEqual(prepared[0]["sample"], "0")
        self.assertEqual(prepared[1]["sample"], "520")
        self.assertEqual(prepared[1]["source_sample"], "500")
        self.assertEqual(prepared[1]["source_onset"], "1")
        self.assertEqual(prepared[1]["source_event_row"], "2")
        self.assertEqual(prepared[1]["value"], "S 10")
        with self.assertRaisesRegex(ValueError, "already prepared"):
            prepare_events(prepared, 500)
        self.assertEqual(prepare_events(rows, 500)[0], prepared)

    def test_native_mne_matching_does_not_reselect_excluded_labels(self):
        import mne
        import numpy as np
        from mne_bids.read import events_file_to_annotation_kwargs
        rows, _ = prepare_events([event(0), event(.3), event(.6)], 500)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.tsv"
            with path.open("w") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            annotations = events_file_to_annotation_kwargs(path, verbose=False)
        raw = mne.io.RawArray(np.zeros((1, 1000)), mne.create_info(["Cz"], 500, "eeg"), verbose=False)
        raw.set_annotations(mne.Annotations(annotations["onset"], annotations["duration"], annotations["description"]))
        events, ids = mne.events_from_annotations(raw, verbose=False)
        selected = mne.event.match_event_names(ids, ["PLAYER_CRASH_WALL"])
        self.assertEqual(selected, ["PLAYER_CRASH_WALL"])
        self.assertNotEqual(ids["PLAYER_CRASH_WALL"], ids[EXCLUDED_PREFIX + "PLAYER_CRASH_WALL"])
        self.assertEqual(sum(events[:, 2] == ids["PLAYER_CRASH_WALL"]), 2)

    def test_union_preserves_actual_per_run_flags_and_eog(self):
        rows = [dict(name=name, type="n/a", units="n/a") for name in ("Cz", "Pz", "VEOG", "HEOG")]
        policy = dict(per_recording_bads={"sub-001_run-01": ["Cz"], "sub-001_run-02": ["Pz"]},
                      subject_union_bads={"001": ["Cz", "Pz"]})
        run1 = prepare_channels(rows, "sub-001_run-01", policy)
        run2 = prepare_channels(rows, "sub-001_run-02", policy)
        self.assertEqual([r["status"] for r in run1], ["bad", "bad", "good", "good"])
        self.assertEqual([r["status"] for r in run1], [r["status"] for r in run2])
        self.assertEqual(run1[1]["pyprep_bad_in_this_run"], "false")
        self.assertEqual(run2[1]["pyprep_bad_in_this_run"], "true")
        self.assertEqual(run1[2]["type"], "EOG")
        self.assertEqual(rows[0]["type"], "n/a")

    def test_preparation_is_independent_idempotent_and_policy_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, destination = root / "source", root / "prepared"
            eeg = source / "sub-001/eeg"
            eeg.mkdir(parents=True)
            proposals = root / "proposals"
            proposals.mkdir()
            policy = dict(random_state=2026, n_samples=1000, proposal_sha256={},
                          per_recording_bads={}, subject_union_bads={"001": ["Cz"]})
            for run in ("01", "02"):
                base = eeg / f"sub-001_task-ContinuousVideoGamePlay_run-{run}"
                Path(str(base) + "_channels.tsv").write_text("name\ttype\tunits\nCz\tn/a\tn/a\nVEOG\tn/a\tn/a\nHEOG\tn/a\tn/a\n")
                Path(str(base) + "_eeg.json").write_text(json.dumps(dict(SamplingFrequency=500)))
                Path(str(base) + "_eeg.set").write_bytes(b"original EEG content")
                Path(str(base) + "_events.tsv").write_text("onset\tduration\tsample\ttrial_type\tvalue\n0\tn/a\t0\tPLAYER_CRASH_WALL\tS 10\n0.3\tn/a\t150\tPLAYER_CRASH_WALL\tS 10\n")
                Path(str(base) + "_events.json").write_text("{}")
                key = f"sub-001_run-{run}"
                bads = ["Cz"] if run == "01" else []
                proposal = proposals / f"{key}.json"
                proposal.write_text(json.dumps(dict(all_bads=bads, settings=dict(random_state=2026, n_samples=1000))))
                policy["proposal_sha256"][key] = hashlib.sha256(proposal.read_bytes()).hexdigest()
                policy["per_recording_bads"][key] = bads
            policy_path = root / "policy.json"
            policy_path.write_text(json.dumps(policy))
            original_bytes = {p: p.read_bytes() for p in source.rglob("*") if p.is_file()}
            kwargs = dict(analysis_ready=True, policy_path=policy_path,
                          proposal_root=proposals, manifest_path=root / "manifest.json")
            prepare(source, destination, **kwargs)
            events = next(destination.rglob("*_events.tsv"))
            first_bytes, first_mtime = events.read_bytes(), events.stat().st_mtime_ns
            self.assertFalse(events.is_symlink())
            self.assertTrue(next(destination.rglob("*_eeg.set")).is_symlink())
            prepare(source, destination, **kwargs)
            self.assertEqual(events.read_bytes(), first_bytes)
            self.assertEqual(events.stat().st_mtime_ns, first_mtime)
            self.assertEqual(original_bytes, {p: p.read_bytes() for p in original_bytes})
            prepared_bytes = {p: p.read_bytes() for p in destination.rglob("*") if p.is_file() and not p.is_symlink()}
            with self.assertRaisesRegex(ValueError, "already analysis-ready"):
                prepare(source, destination, manifest_path=root / "default_manifest.json")
            self.assertEqual(prepared_bytes, {p: p.read_bytes() for p in prepared_bytes})
            self.assertEqual(len(list(destination.rglob("*_channels.json"))), 2)
            proposal.write_text("{}")
            with self.assertRaisesRegex(ValueError, "detector hash"):
                prepare(source, destination, **kwargs)
            self.assertEqual(events.read_bytes(), first_bytes)


if __name__ == "__main__":
    unittest.main()
