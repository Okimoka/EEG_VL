"""Fixed coordinate repair for this dataset's shared BESA spherical template.

Only the coordinate frame changes. Landmarks are from the matching template,
not participant digitization. No template download is needed during preparation.
"""
import csv
import hashlib
import io

import mne
import numpy as np

BESA_URL = "https://raw.githubusercontent.com/sccn/dipfit/master/standard_BESA/standard-10-5-cap385.elp"
BESA_SHA256 = "b94ef2b246b61ff60bb5c29372ce5b7ee3e27074a7cc542920f27afa46729697"
# All 34 originals match the named template to saved precision. Refuse an
# unaudited geometry or an accidental second transformation.
SOURCE_ELECTRODES_SHA256 = "d7aea104c5cfe494bf22a2a3f1e845fecd10d65d1999641db4c77260db91786b"
BESA_LANDMARKS = {"LPA": (-120.03, 0.0, 85.0),
                  "RPA": (120.03, 0.0, 85.0),
                  "NAS": (114.03, 90.0, 85.0)}
ALS_TO_RAS = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])


def template_landmarks_ras():
    """Convert the three BESA spherical landmark rows to RAS metres."""
    positions = {}
    for name, (azimuth, horizontal, radius_mm) in BESA_LANDMARKS.items():
        az = np.deg2rad(horizontal if azimuth >= 0 else 180 + horizontal)
        polar = np.deg2rad(abs(azimuth))
        positions[name] = radius_mm * 1e-3 * np.array(
            [np.sin(polar) * np.cos(az), np.sin(polar) * np.sin(az), np.cos(polar)])
    return positions


def ras_to_head_transform():
    landmarks = template_landmarks_ras()
    montage = mne.channels.make_dig_montage(
        nasion=landmarks["NAS"], lpa=landmarks["LPA"], rpa=landmarks["RPA"],
        coord_frame="unknown")
    return mne.channels.compute_native_head_t(montage)["trans"]



def topomap_sphere():
    """The original 85 mm template sphere, expressed in the output head frame."""
    return np.r_[ras_to_head_transform()[:3, 3], .085]


def corrected_electrodes_tsv(path):
    """Return a table in the MNE/CapTrak head frame, in metres."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_ELECTRODES_SHA256:
        raise ValueError(f"Electrode geometry differs from the audited BESA source: {path}")
    rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))
    native = np.array([[float(row[c]) for c in ("x", "y", "z")] for row in rows])
    ras = native @ ALS_TO_RAS.T * 1e-3
    head = mne.transforms.apply_trans(ras_to_head_transform(), ras)
    for row, position in zip(rows, head):
        for axis, value in zip(("x", "y", "z"), position):
            row[axis] = format(value, ".12g")
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return text.getvalue()


def corrected_coordsystem(original):
    if original.get("EEGCoordinateUnits") != "mm":
        raise ValueError("Expected the audited source electrode coordinates in millimetres")
    result = dict(original)
    transform = ras_to_head_transform()
    landmarks = {name: mne.transforms.apply_trans(transform, pos).tolist()
                 for name, pos in template_landmarks_ras().items()}
    result.update(
        EEGCoordinateSystem="CapTrak", EEGCoordinateUnits="m",
        EEGCoordinateSystemDescription=(
            "RAS head coordinates: +X right, +Y anterior, +Z superior; origin at "
            "the midpoint of the BESA template LPA/RPA, nasion defining +Y. "
            "Supplied spherical BESA positions were rigidly transformed using "
            "the matching standard-10-5-cap385.elp template landmarks. CapTrak "
            "names the coordinate convention only; no CapTrak device or "
            "individual digitization was used. All participants share template "
            "positions. VEOG/HEOG rows remain transformed source placeholders, "
            "not scalp EEG electrodes."),
        AnatomicalLandmarkCoordinates=landmarks,
        AnatomicalLandmarkCoordinateSystem="CapTrak",
        AnatomicalLandmarkCoordinateUnits="m",
        AnatomicalLandmarkCoordinateSystemDescription=(
            "Inferred BESA template landmarks in the same head frame as the "
            "electrodes; not measured participant anatomy."))
    return result


def montage_provenance():
    return dict(
        method="BESA template-fiducial rigid transform; no electrode substitution or scaling",
        source_electrodes_sha256=SOURCE_ELECTRODES_SHA256,
        template_url=BESA_URL, template_sha256=BESA_SHA256,
        template_landmarks_spherical=BESA_LANDMARKS,
        spherical_columns=["signed_polar_degrees", "horizontal_degrees", "radius_mm"],
        source_frame="ALS spherical-template origin", source_units="mm",
        output_frame="CapTrak / MNE head", output_units="m",
        source_ALS_to_RAS_rotation=ALS_TO_RAS.tolist(),
        RAS_metres_to_head_metres=ras_to_head_transform().tolist(),
        landmarks_are_measured=False)
