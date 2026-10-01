#!/usr/bin/env python3
"""Desktop review of fixed PyPREP proposals, with decisions in channel_review.xlsx.
Fully written by an LLM
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
REVIEW_FOLDER = Path(__file__).resolve().parent.name
KEY = re.compile(r"sub-([A-Za-z0-9]+)_run-([A-Za-z0-9]+)\Z")
LABELS = {
    "bad_by_psd": "PSD",
    "bad_by_ransac": "RANSAC",
    "bad_by_correlation": "correlation",
    "bad_by_deviation": "deviation",
    "bad_by_hf_noise": "high-frequency noise",
    "bad_by_dropout": "dropout",
    "bad_by_flat": "flat signal",
    "bad_by_nan": "NaN",
    "bad_by_SNR": "SNR",
    "bad_by_manual": "prior manual flag",
}


def flag_reasons(proposal, channel):
    return [LABELS.get(key, key.removeprefix("bad_by_"))
            for key, names in proposal.get("bads", {}).items()
            if key != "bad_all" and channel in names]


def flagged_channels(proposal):
    return set(proposal.get("all_bads", [])) | {
        channel for names in proposal.get("bads", {}).values() for channel in names
    }


@dataclass(frozen=True)
class Bookmark:
    label: str
    start: float
    detail: str = ""


def bookmarks(proposal, channel, duration, diagnostics_path=None):
    """Rank actual failing windows; never pretend a global PSD flag has a trigger."""
    flagged_criteria = {key for key, names in proposal.get("bads", {}).items()
                        if key != "bad_all" and channel in names}
    groups = {}
    for window in proposal.get("windows", []):
        if window.get("channel") != channel:
            continue
        criterion = window["criterion"]
        if criterion not in {"bad_by_correlation", "bad_by_ransac", "bad_by_dropout"}:
            continue
        # These exported windows represent threshold failures, even if too few
        # failed for the channel to receive a whole-record flag.
        groups.setdefault(criterion, []).append(window)
    result = []
    for criterion in sorted(groups, key=lambda key: (key not in flagged_criteria, key)):
        windows = groups[criterion]
        descending = criterion == "bad_by_dropout"
        def score(window):
            value = window.get("value")
            if value is None:
                return float("inf")
            return -float(value) if descending else float(value)
        chosen = []
        for window in sorted(windows, key=score):
            start = float(window["start"])
            if any(abs(start - other) < 20 for other in chosen):
                continue
            chosen.append(start)
            value = window.get("value")
            value_text = f"; score {value:.3g}" if isinstance(value, (float, int)) else ""
            label = (f"{LABELS.get(criterion, criterion)} failing window: "
                     f"{start:g}–{float(window['stop']):g} s{value_text}")
            detail = window.get("detail", "")
            if criterion not in flagged_criteria:
                detail = "This channel was not globally flagged by this detector. " + detail
            result.append(Bookmark(label, max(0.0, start - 5), detail))
            if len(chosen) == 3:
                break
    contexts = []
    if diagnostics_path is not None and Path(diagnostics_path).exists():
        import numpy as np
        # These are navigation aids for global tests, not their actual triggers.
        fields = []
        if flagged_criteria & {"bad_by_deviation", "bad_by_psd"}:
            fields.append(("bad_by_deviation__channel_amplitudes", "High-amplitude", 1e6, "µV"))
        if "bad_by_hf_noise" in flagged_criteria:
            fields.append(("bad_by_hf_noise__noise_levels", "High-frequency noise", 1.0, "ratio"))
        names = proposal.get("eeg_ch_names", [])
        if channel in names:
            column = names.index(channel)
            with np.load(diagnostics_path, allow_pickle=False) as diagnostics:
                for field, label, scale, unit in fields:
                    if field not in diagnostics:
                        continue
                    values = diagnostics[field][:, column]
                    chosen = []
                    for index in np.argsort(-values):
                        if not np.isfinite(values[index]) or any(abs(int(index) - old) < 20 for old in chosen):
                            continue
                        chosen.append(int(index))
                        contexts.append(Bookmark(
                            f"{label} context (not a trigger): {index}–{index + 1} s "
                            f"({values[index] * scale:.3g} {unit})",
                            max(0.0, float(index) - 5),
                            "A high-value one-second diagnostic interval, selected only for inspection. "
                            "The channel flag uses a whole-record statistic; this interval did not "
                            "independently trigger it. In particular, amplitude is not the PSD test.",
                        ))
                        if len(chosen) == 3:
                            break
    comparison = [Bookmark("Overview: start of recording", 0.0)]
    comparison.extend(Bookmark(f"Comparison: {int(fraction * 100)}% of recording",
                               max(0.0, duration * fraction - 5))
                      for fraction in (0.1, 0.5, 0.9))
    # Prioritize contextual extremes for global flags, then neutral comparisons.
    localized_flag = bool(flagged_criteria & groups.keys())
    return result + contexts + comparison if localized_flag else contexts + comparison + result


def spatial_order(raw, channel):
    """Selected EEG channel, then every other EEG by Euclidean montage distance."""
    import numpy as np
    eeg = [i for i, kind in enumerate(raw.get_channel_types()) if kind == "eeg"]
    selected = raw.ch_names.index(channel)
    origin = raw.info["chs"][selected]["loc"][:3]
    def distance(index):
        position = raw.info["chs"][index]["loc"][:3]
        if not np.all(np.isfinite(origin)) or not np.all(np.isfinite(position)):
            return float("inf")
        return float(np.linalg.norm(position - origin))
    return [selected] + sorted((index for index in eeg if index != selected), key=distance)


def configure_runtime(root):
    runtime = root / "logs" / "review" / "runtime"
    for key, suffix in (("MPLCONFIGDIR", "mpl"), ("_MNE_FAKE_HOME_DIR", "mne")):
        folder = runtime / suffix
        folder.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault(key, str(folder))
    os.environ.setdefault("QT_API", "pyqt6")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def reviewer_class():
    """Defer GUI imports so --help works without a display or scientific imports."""
    import mne
    import numpy as np
    from PyQt6 import QtCore, QtWidgets

    class Reviewer(QtWidgets.QMainWindow):
        def __init__(self, store, root=ROOT, key=None):
            super().__init__()
            self.store, self.root = store, Path(root)
            self.raw = self.browser = self.eog_browser = None
            self.current_ch = self.key = None
            self.loading = False
            self.choices = {}
            self.baseline = {}
            self.pending_overrides, self.pending_notes = {}, {}
            self.loaded_mode = 0
            self.setWindowTitle("EEG channel review")
            self.resize(620, 780)
            panel = QtWidgets.QWidget()
            self.setCentralWidget(panel)
            layout = QtWidgets.QVBoxLayout(panel)

            self.recording = QtWidgets.QComboBox()
            for record_key, proposal in sorted(store.recordings.items()):
                subject, run = KEY.fullmatch(record_key).groups()
                task = "OddGamble" if int(run) == 1 else "Gameplay"
                self.recording.addItem(f"{subject} — {task}", record_key)
            if key is not None:
                index = self.recording.findData(key)
                if index < 0:
                    raise ValueError(f"Recording is not available: {key}")
                self.recording.setCurrentIndex(index)
            layout.addWidget(self.recording)

            self.mode = QtWidgets.QComboBox()
            self.mode.addItems(["0.1–100 Hz + 60 Hz notch (EEG)", "Original signal"])
            layout.addWidget(self.mode)
            self.scope = QtWidgets.QComboBox()
            self.scope.addItems(["Flagged EEG channels", "All EEG channels"])
            layout.addWidget(self.scope)
            self.channels = QtWidgets.QListWidget()
            layout.addWidget(self.channels, 2)
            self.reasons = QtWidgets.QLabel()
            self.reasons.setWordWrap(True)
            layout.addWidget(self.reasons)
            self.windows = QtWidgets.QComboBox()
            self.windows.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            self.windows.setMinimumContentsLength(32)
            layout.addWidget(self.windows)

            row = QtWidgets.QHBoxLayout()
            layout.addLayout(row)
            self.open_button = QtWidgets.QPushButton("Open / jump in MNE")
            self.next_button = QtWidgets.QPushButton("Next channel")
            self.eog_button = QtWidgets.QPushButton("EOG comparison")
            for button in (self.open_button, self.next_button, self.eog_button):
                row.addWidget(button)
            self.override = QtWidgets.QCheckBox("Manual override")
            self.override.setToolTip("Checked: invert the PyPREP decision. Unchecked: accept it.")
            self.final_decision = QtWidgets.QLabel()
            row = QtWidgets.QHBoxLayout()
            row.addWidget(self.override)
            row.addWidget(self.final_decision)
            layout.addLayout(row)
            self.note = QtWidgets.QPlainTextEdit()
            self.note.setPlaceholderText("Optional decision note and timestamps")
            self.note.setMaximumHeight(110)
            layout.addWidget(self.note)
            row = QtWidgets.QHBoxLayout()
            layout.addLayout(row)
            self.save_button = QtWidgets.QPushButton("Save draft")
            self.export_button = QtWidgets.QPushButton("Export this interval")
            self.reload_button = QtWidgets.QPushButton("Reload spreadsheet")
            for button in (self.save_button, self.export_button, self.reload_button):
                row.addWidget(button)
            row = QtWidgets.QHBoxLayout()
            layout.addLayout(row)
            self.reviewed_button = QtWidgets.QPushButton("Mark recording reviewed")
            self.review_status = QtWidgets.QLabel()
            row.addWidget(self.reviewed_button)
            row.addWidget(self.review_status)
            self.message = QtWidgets.QLabel()
            self.message.setWordWrap(True)
            layout.addWidget(self.message)

            self.recording.currentIndexChanged.connect(lambda: self.guarded(self.load_recording))
            self.mode.currentIndexChanged.connect(lambda: self.guarded(self.load_mode))
            self.scope.currentIndexChanged.connect(lambda: self.guarded(self.populate_channels))
            self.channels.currentRowChanged.connect(lambda: self.guarded(self.select_channel))
            self.override.toggled.connect(lambda value: self.guarded(lambda: self.change_override(value)))
            self.note.textChanged.connect(self.note_changed)
            self.open_button.clicked.connect(lambda: self.guarded(self.open_browser))
            self.eog_button.clicked.connect(lambda: self.guarded(self.open_eog))
            self.next_button.clicked.connect(self.next_channel)
            self.save_button.clicked.connect(lambda: self.guarded(self.save))
            self.export_button.clicked.connect(lambda: self.guarded(self.export_evidence))
            self.reload_button.clicked.connect(lambda: self.guarded(self.reload_workbook))
            self.reviewed_button.clicked.connect(lambda: self.guarded(self.mark_reviewed))
            self.load_recording()

        def guarded(self, function):
            try:
                function()
            except Exception as exc:
                self.loading = False
                QtWidgets.QMessageBox.critical(self, "Channel review", str(exc))

        @property
        def proposal(self):
            return self.store.recordings[self.key]

        def note_changed(self):
            if self.loading or not self.current_ch:
                return
            note = self.note.toPlainText()
            if note != self.choices[self.current_ch].get("note", ""):
                self.pending_notes[self.current_ch] = note
                self.choices[self.current_ch]["note"] = note
                self.review_status.setText("Draft")

        def update_item(self, channel):
            reason = ", ".join(flag_reasons(self.proposal, channel))
            text = f"{channel} — flagged by {reason}" if reason else f"{channel} — not flagged"
            if self.choices[channel]["override"]:
                text += " · override"
            for row in range(self.channels.count()):
                item = self.channels.item(row)
                if item.data(QtCore.Qt.ItemDataRole.UserRole) == channel:
                    item.setText(text)
                    break

        def sync_current_widgets(self):
            if not self.current_ch:
                return
            choice = self.choices[self.current_ch]
            self.override.blockSignals(True)
            self.override.setChecked(choice["override"])
            self.override.blockSignals(False)
            self.final_decision.setText("Final decision: " + ("mark bad" if choice["final_bad"] else "keep"))

        def capture_browser(self):
            if self.raw is None or self.browser is None:
                return
            bads = set(self.browser.mne.info["bads"])
            flags = flagged_channels(self.proposal)
            for channel, choice in self.choices.items():
                final_bad = channel in bads
                override = final_bad != (channel in flags)
                if override != choice["override"]:
                    choice.update(override=override, final_bad=final_bad)
                    self.pending_overrides[channel] = override
                    self.update_item(channel)
                    self.review_status.setText("Draft")
            self.raw.info["bads"] = [ch for ch, value in self.choices.items() if value["final_bad"]]
            self.sync_current_widgets()

        def save(self):
            if self.key is None:
                return
            self.capture_browser()
            if self.pending_overrides or self.pending_notes:
                fresh = self.store.choices(self.key)
                conflicts = []
                for field, changes in (("override", self.pending_overrides), ("note", self.pending_notes)):
                    for channel, value in changes.items():
                        if fresh[channel][field] not in (self.baseline[channel][field], value):
                            conflicts.append(f"{channel} {field}")
                if conflicts:
                    raise ValueError("The spreadsheet was edited outside the application: "
                                     + ", ".join(conflicts)
                                     + ". Use Reload spreadsheet to read those changes; "
                                     "your pending edits have not been written.")
                self.store.save(self.key, dict(self.pending_overrides), dict(self.pending_notes), reviewed=False)
                for field, changes in (("override", self.pending_overrides), ("note", self.pending_notes)):
                    for channel, value in changes.items():
                        self.baseline[channel][field] = value
                self.pending_overrides.clear()
                self.pending_notes.clear()
            self.review_status.setText("Reviewed" if self.store.is_reviewed(self.key) else "Draft")
            self.message.setText(f"Saved draft: {self.store.path.name} · {self.key}")

        def _read_raw(self):
            if self.mode.currentIndex() == 0:
                path = self.root / self.proposal["input"]
                raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
            else:
                from mne_bids import BIDSPath, read_raw_bids
                subject, run = KEY.fullmatch(self.key).groups()
                path = BIDSPath(root=self.root / "prepared_bids", subject=subject,
                                run=run, session=self.proposal.get("session"),
                                task=self.proposal.get("task", "ContinuousVideoGamePlay"),
                                datatype="eeg", suffix="eeg", extension=".set")
                raw = read_raw_bids(path, verbose="ERROR")
            # Hide dense task-event labels in this channel-quality view only.
            raw.set_annotations(mne.Annotations([], [], []))
            raw.info["bads"] = [ch for ch, value in self.choices.items() if value["final_bad"]]
            return raw

        def close_browsers(self):
            self.capture_browser()
            for name in ("browser", "eog_browser"):
                browser = getattr(self, name)
                if browser is not None:
                    setattr(self, name, None)
                    browser.close()

        def load_recording(self):
            if self.loading:
                return
            requested = self.recording.currentData()
            try:
                self.save()
            except Exception:
                self.recording.blockSignals(True)
                self.recording.setCurrentIndex(self.recording.findData(self.key))
                self.recording.blockSignals(False)
                raise
            self.close_browsers()
            if self.raw is not None:
                self.raw.close()
            self.key = requested
            self.current_ch = None
            self.choices = self.store.choices(self.key)
            self.baseline = deepcopy(self.choices)
            self.raw = self._read_raw()
            self.loaded_mode = self.mode.currentIndex()
            self.populate_channels()
            self.review_status.setText("Reviewed" if self.store.is_reviewed(self.key) else "Draft")
            self.message.clear()

        def load_mode(self):
            if self.loading or self.key is None:
                return
            reopen = self.browser is not None
            try:
                self.save()
                new_raw = self._read_raw()
            except Exception:
                self.mode.blockSignals(True)
                self.mode.setCurrentIndex(self.loaded_mode)
                self.mode.blockSignals(False)
                raise
            self.close_browsers()
            self.raw.close()
            self.raw = new_raw
            self.loaded_mode = self.mode.currentIndex()
            if reopen:
                self.open_browser()

        def populate_channels(self):
            if self.loading or self.raw is None:
                return
            self.save()
            previous = self.current_ch
            self.loading = True
            flags = flagged_channels(self.proposal)
            names = [name for name in self.raw.ch_names if name in self.choices]
            if self.scope.currentIndex() == 0:
                names = [name for name in names if name in flags]
            if not names:
                self.scope.setCurrentIndex(1)
                names = [name for name in self.raw.ch_names if name in self.choices]
            names.sort(key=lambda name: name not in flags)
            self.channels.clear()
            for name in names:
                item = QtWidgets.QListWidgetItem()
                item.setData(QtCore.Qt.ItemDataRole.UserRole, name)
                self.channels.addItem(item)
                self.update_item(name)
            self.current_ch = None
            self.loading = False
            self.channels.setCurrentRow(names.index(previous) if previous in names else 0)

        def select_channel(self):
            if self.loading:
                return
            self.save()
            item = self.channels.currentItem()
            if item is None:
                return
            self.current_ch = item.data(QtCore.Qt.ItemDataRole.UserRole)
            self.loading = True
            self.note.setPlainText(self.choices[self.current_ch].get("note", ""))
            self.sync_current_widgets()
            reason = ", ".join(flag_reasons(self.proposal, self.current_ch))
            self.reasons.setText(f"{self.current_ch} — flagged by {reason}" if reason
                                 else f"{self.current_ch} — not flagged by PyPREP")
            self.windows.clear()
            diagnostics = self.root / "logs" / "pyprep" / f"{self.key}_diagnostics.npz"
            for bookmark in bookmarks(self.proposal, self.current_ch,
                                      self.raw.n_times / self.raw.info["sfreq"], diagnostics):
                self.windows.addItem(bookmark.label, bookmark.start)
                self.windows.setItemData(self.windows.count() - 1, bookmark.detail,
                                         QtCore.Qt.ItemDataRole.ToolTipRole)
            self.loading = False
            if self.browser is not None:
                self.open_browser()

        def change_override(self, value):
            if self.loading or self.current_ch is None:
                return
            self.capture_browser()
            final_bad = (self.current_ch in flagged_channels(self.proposal)) != value
            self.choices[self.current_ch].update(override=value, final_bad=final_bad)
            self.pending_overrides[self.current_ch] = value
            bads = [ch for ch, choice in self.choices.items() if choice["final_bad"]]
            self.raw.info["bads"] = bads
            if self.browser is not None:
                self.browser.mne.info["bads"] = list(bads)
                self.browser._redraw()
            self.sync_current_widgets()
            self.update_item(self.current_ch)
            self.save()

        def next_channel(self):
            self.channels.setCurrentRow(min(self.channels.currentRow() + 1, self.channels.count() - 1))

        def open_browser(self):
            if not self.current_ch:
                return
            self.close_browsers()
            duration = min(15.0, self.raw.n_times / self.raw.info["sfreq"])
            start = min(float(self.windows.currentData() or 0), max(0.0, self.raw.times[-1] - duration))
            order = spatial_order(self.raw, self.current_ch)
            mne.viz.set_browser_backend("qt")
            self.browser = self.raw.plot(
                order=order, n_channels=min(12, len(order)), duration=duration,
                start=start, scalings=dict(eeg=50e-6, eog=5e-3), bad_color="#d65f5f",
                title=f"{self.recording.currentText()} · {self.current_ch} · {self.mode.currentText()}",
                show=False, block=False, precompute=False, show_options=False,
                proj=False, remove_dc=True, group_by="original", decim=1,
            )
            self.browser.installEventFilter(self)
            self.browser.resize(1250, 780)
            self.browser.show()

        def open_eog(self):
            if self.eog_browser is not None:
                self.eog_browser.close()
            eeg_context = [name for name in ("Fp1", "Fp2") if name in self.raw.ch_names]
            names = [name for name, kind in zip(self.raw.ch_names, self.raw.get_channel_types())
                     if kind == "eog"] + eeg_context
            order = [self.raw.ch_names.index(name) for name in names]
            if not order:
                return
            start = float(self.browser.mne.t_start) if self.browser is not None else float(self.windows.currentData() or 0)
            mne.viz.set_browser_backend("qt")
            # A separate copy keeps EOG-window markings out of the EEG decisions.
            self.eog_raw = self.raw.copy()
            self.eog_browser = self.eog_raw.plot(
                order=order, n_channels=len(order), duration=15, start=start,
                scalings=dict(eeg=50e-6, eog=5e-3), show=False, block=False,
                precompute=False, show_options=False, proj=False, remove_dc=True,
                group_by="original", decim=1, title=f"{self.recording.currentText()} · EOG comparison",
            )
            self.eog_browser.installEventFilter(self)
            self.eog_browser.resize(1150, 550)
            self.eog_browser.show()

        def eventFilter(self, obj, event):
            if event.type() == QtCore.QEvent.Type.Close:
                if obj is self.browser:
                    self.capture_browser()
                    self.browser = None
                    self.guarded(self.save)
                elif obj is self.eog_browser:
                    self.eog_browser = None
            return super().eventFilter(obj, event)

        def reload_workbook(self):
            self.capture_browser()
            if self.pending_overrides or self.pending_notes:
                response = QtWidgets.QMessageBox.question(
                    self, "Reload spreadsheet", "Discard pending application edits and reload the spreadsheet?",
                    QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
                    QtWidgets.QMessageBox.StandardButton.No,
                )
                if response != QtWidgets.QMessageBox.StandardButton.Yes:
                    return
            self.close_browsers()
            self.store.reload()
            self.choices = self.store.choices(self.key)
            self.baseline = deepcopy(self.choices)
            self.pending_overrides.clear()
            self.pending_notes.clear()
            self.raw.info["bads"] = [ch for ch, choice in self.choices.items() if choice["final_bad"]]
            self.populate_channels()
            self.review_status.setText("Reviewed" if self.store.is_reviewed(self.key) else "Draft")

        def mark_reviewed(self):
            self.save()
            self.store.save(self.key, {}, reviewed=True)
            self.review_status.setText("Reviewed")
            self.message.setText(f"Saved draft: {self.store.path.name} · {self.key} · reviewed")

        def export_evidence(self):
            if self.browser is None:
                raise ValueError("Open the MNE browser at the interval to export first.")
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            self.save()
            start, duration = float(self.browser.mne.t_start), float(self.browser.mne.duration)
            picks = list(self.browser.mne.picks)
            sfreq = self.raw.info["sfreq"]
            first, last = max(0, int(round(start * sfreq))), min(self.raw.n_times, int(round((start + duration) * sfreq)))
            data = self.raw.get_data(picks=picks, start=first, stop=last) * 1e6
            centered = data - np.mean(data, axis=1, keepdims=True)
            limit = max(50.0, float(np.max(np.abs(centered))) * 1.05)
            fig = Figure(figsize=(12, max(4, len(picks) * 0.55)))
            FigureCanvasAgg(fig)
            axes = fig.subplots(len(picks), 1, sharex=True, squeeze=False).ravel()
            for ax, signal, index in zip(axes, centered, picks):
                channel = self.raw.ch_names[index]
                ax.plot(np.arange(first, last) / sfreq, signal, lw=0.6,
                        color="#aa3030" if channel == self.current_ch else "#265c80")
                ax.set_ylim(-limit, limit)
                ax.set_ylabel(channel, rotation=0, labelpad=24, fontsize=9)
            axes[-1].set_xlabel("Recording time (s)")
            fig.suptitle(f"{self.recording.currentText()} · {self.current_ch} · {self.mode.currentText()}\n"
                         f"Acquisition reference; interval means removed; common scale ±{limit:.0f} µV", fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            folder = self.root / "figures" / "channel_review"
            folder.mkdir(parents=True, exist_ok=True)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            path = folder / f"{self.key}_{self.current_ch}_{stamp}.png"
            fig.savefig(path, dpi=160)
            metadata = dict(recording=self.key, channel=self.current_ch,
                            start_s=first / sfreq, stop_s=last / sfreq,
                            displayed_channels=[self.raw.ch_names[index] for index in picks],
                            signal_view=self.mode.currentText(), reference="acquisition",
                            centering="interval mean removed for display only", scale_uV=limit,
                            choice=self.choices[self.current_ch], proposal_settings=self.proposal.get("settings", {}),
                            saved_utc=datetime.now(timezone.utc).isoformat())
            path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
            self.message.setText(f"Saved draft: {self.store.path.name} · figure {path.name}")

        def closeEvent(self, event):
            try:
                self.save()
                self.close_browsers()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Unsaved channel review", str(exc))
                event.ignore()
                return
            if self.raw is not None:
                self.raw.close()
            event.accept()

    return Reviewer


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", default="001")
    parser.add_argument("--run", choices=("01", "02"), default="01")
    parser.add_argument("--root", type=Path, default=ROOT, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    configure_runtime(root)
    import fcntl
    import mne
    from PyQt6 import QtWidgets
    from review_store import ReviewStore
    mne.set_log_level("ERROR")
    lock_path = root / REVIEW_FOLDER / ".review.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.exit(1, "A channel-review application is already open.\n")
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
        try:
            store = ReviewStore(root=root)
            if not store.path.exists():
                store.initialize()
            window = reviewer_class()(store, root=root, key=f"sub-{args.subject.zfill(3)}_run-{args.run}")
        except Exception as exc:
            QtWidgets.QMessageBox.critical(None, "Cannot open channel review", str(exc))
            return 1
        window.show()
        window.guarded(window.open_browser)
        return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
