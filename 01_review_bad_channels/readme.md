# Channel inspection

Run after preliminary filtering and PyPREP (see the root readme):

```bash
python 01_review_bad_channels/review_store.py
./01_review_bad_channels/review-channels --subject 001 --run 01
```

Select a recording and a channel. The browser shows that electrode together with nearby electrodes, chosen by distance in the corrected montage; FT7 is compared with T7 because T7 is nearby, not because PyPREP tested only that pair. Switch between original and filtered EEG and use the interval selector for detector-specific evidence. Correlation/dropout use one-second windows; RANSAC uses five-second windows. Amplitude/PSD flags may refer to the whole recording. MNE may reduce display resolution when many samples fit on screen; the tool does not average neighbouring channels into one signal.

The JSONs contain failing-window bookmarks. Fresh PyPREP runs also save detailed diagnostic NPZ arrays; these large arrays are not included in the submission. The local workbook records observations and optional overrides, but the final project policy explicitly accepts all PyPREP flags. `scripts/channel_policy.py` exports that policy independently of exploratory workbook edits.
