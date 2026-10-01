Run from the submission root:

```sh
python analysis/prepare.py
python analysis/train.py
python analysis/results.py
```

Inputs: the existing cleaned epochs and continuous recordings in `../artifacts/postica/`.
Generated EEG, models and numeric results go to `../artifacts/analysis_final/`.
Only the report's four tables and transfer figure are written into this submission.
Training uses an NVIDIA GPU. Settings are written directly in the scripts.
