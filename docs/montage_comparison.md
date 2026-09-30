# Geometry used in preparation

All 34 supplied `*_electrodes.tsv` files have SHA256 `d7aea104c5cfe494bf22a2a3f1e845fecd10d65d1999641db4c77260db91786b`. They match the BESA spherical `standard-10-5-cap385.elp` template named in the supplied MATLAB scripts, to the files' 0.01 mm rounding. They are not the geometry of MNE's `standard_1005` or `standard_1020`.

`scripts/_montage.py` retains every coordinate and applies an ALS-to-RAS rotation followed by a rigid transform based on the matching template nasion/LPA/RPA. Coordinates become metres in MNE head space. The BIDS sidecar calls this convention `CapTrak`; no CapTrak recording or individual digitization is implied. Preparation refuses a different input-table hash, preventing accidental double transformation or substitution of an unaudited montage.

The template's URL, hash, three spherical landmark rows and the transformation matrix are recorded by `montage_provenance()` and written into each preparation receipt. The needed landmarks are embedded in the script, so preparation needs no external template download. EEG/EOG counts are 63/2, with CPz as acquisition reference but no separately recorded CPz channel. EOG coordinate rows remain placeholders and do not enter scalp-channel detection or ICA.
