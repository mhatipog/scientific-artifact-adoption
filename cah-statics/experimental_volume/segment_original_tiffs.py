"""Independent phase segmentation of Yan et al. (2025) original Q1 water-air TIFFs.

Source: https://doi.org/10.18419/DARUS-4683 (CC BY 4.0).
Run manually, one ~200/345 MB download at a time. Never commit the raw TIFF.
Classification is a preliminary 2D volume proxy; bright unclassified pixels
include corner gas, thin films and optical artifacts. No pressure clock is inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter

PITCH_UM = 2.4859
PORE_WIDTH_UM = 800.0
DEPTH_UM = 120.0
PORE_VOLUME_UL = PORE_WIDTH_UM**2 * DEPTH_UM / 1e9
PUMP_UL_S = 0.01 * 1000 / 3600
MAX_BYTES = 430_000_000
FILES = {
    60: dict(id=351428, md5="3dbfc99446d61b89c905057cd360d5b7",
             frames=381, dt=0.03332185745239258),
    120: dict(id=351438, md5="1beeb40ae47b8c7299aba034e40c5c5f",
              frames=658, dt=0.033332135528326035),
}
SETTINGS = ((100, .4), (180, .4), (180, .5), (180, .6), (260, .6))
BASE = "https://darus.uni-stuttgart.de/api/access/datafile/"


def crossing(time, signal, value):
    """First crossing of monotone upper envelope, with linear interpolation."""
    env = np.maximum.accumulate(np.asarray(signal, dtype=float))
    ix = np.flatnonzero(env >= value)
    if len(ix) == 0:
        return None
    j = int(ix[0])
    if j == 0:
        return float(time[0])
    rise = env[j] - env[j - 1]
    if rise <= 0:
        raise AssertionError("Nonpositive rise at first crossing")
    return float(time[j - 1] + (value - env[j - 1])
                 * (time[j] - time[j - 1]) / rise)


def find_pore(first):
    """Locate opposing dark wall lines at the independently specified pore span.

    Search windows apply ONLY to the two published 512x512 TIFF stacks.
    In particular this is not a generic automatic microchannel detector.
    """
    pixels = round(PORE_WIDTH_UM / PITCH_UM)
    s = gaussian_filter(first, 1.2)
    row = np.median(s[:, 180:335], axis=1)
    yc = np.arange(20, 110)
    y0 = int(yc[np.argmin(row[yc] + row[yc + pixels])])
    y1 = y0 + pixels
    col = np.median(s[y0 + 70:y1 - 70, :], axis=0)
    xc = np.arange(60, 150)
    x0 = int(xc[np.argmin(col[xc] + col[xc + pixels])])
    x1 = x0 + pixels
    if not (90 <= x0 <= 105 and 25 <= y0 <= 70
            and abs((x1 - x0) * PITCH_UM - PORE_WIDTH_UM) < 3):
        raise RuntimeError(f"Pore-wall detector outside predeclared window: {x0,y0,x1,y1}")
    return x0, y0, x1, y1


def read_original(angle, local_tiff=None, allow_download=False):
    spec = FILES[angle]
    if local_tiff:
        raw = Path(local_tiff).read_bytes()
    else:
        if not allow_download:
            raise ValueError("Explicit --download-original required; raw source is hundreds of MB")
        request = Request(BASE + str(spec["id"]), headers={
            "User-Agent": "IndependentCAHScientificAudit/1.0", "Accept": "image/tiff,*/*"})
        with urlopen(request, timeout=180) as f:
            raw = f.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError("TIFF exceeded 430 MB cap")
    md5 = hashlib.md5(raw).hexdigest()
    if md5 != spec["md5"]:
        raise ValueError(f"Original DaRUS MD5 mismatch ({md5})")
    return raw


def reconstruct(angle, raw):
    spec = FILES[angle]
    image = Image.open(io.BytesIO(raw))
    if image.size != (512, 512) or image.n_frames != spec["frames"]:
        raise ValueError("Unexpected original TIFF dimensions/frame count")
    description = str(image.tag_v2.get(270, ""))
    m = re.search(r"finterval\s*=\s*([0-9.]+)", description)
    dt = float(m.group(1)) if m else spec["dt"]
    if abs(dt - spec["dt"]) > 1e-5:
        raise ValueError("TIFF frame interval disagrees with recorded metadata")

    def frame(i):
        image.seek(int(i))
        return np.asarray(image, dtype=np.float32).copy()

    first, last = frame(0), frame(image.n_frames - 1)
    x0, y0, x1, y1 = find_pore(first)
    roi = np.zeros(first.shape, dtype=bool)
    # Exclude a predeclared three-pixel wall strip; this can miss corner films.
    roi[y0 + 3:y1 - 3, x0 + 3:x1 - 3] = True
    denominator = int(roi.sum())
    a, z = gaussian_filter(first, 1.0), gaussian_filter(last, 1.0)
    contrast = a - z
    exterior = np.zeros(first.shape, dtype=bool)
    exterior[:max(0, y0 - 10), 180:330] = True
    exterior[min(512, y1 + 10):, 180:330] = True
    if not exterior.any():
        raise RuntimeError("No static exterior available for illumination correction")
    active = {(threshold, fraction): roi & (contrast > threshold)
              for threshold, fraction in SETTINGS}
    curves = {key: np.zeros(image.n_frames) for key in SETTINGS}
    for j in range(image.n_frames):
        current = gaussian_filter(frame(j), 1.0)
        drift = float(np.median((current - a)[exterior]))
        darkening = a - current + drift
        for key, eligible in active.items():
            curves[key][j] = np.count_nonzero(
                eligible & (darkening > key[1] * contrast)) / denominator

    t = np.arange(image.n_frames) * dt
    nominal_final = 1.0 if angle == 60 else 0.92
    results = {}
    for key, sat in curves.items():
        t10 = crossing(t, sat, .1 * nominal_final)
        t90 = crossing(t, sat, .9 * nominal_final)
        results[f"{key[0]}_{key[1]:.1f}"] = {
            "final_segmentation_fraction": round(float(sat[-1]), 5),
            "t10_clip_s": t10, "t90_clip_s": t90,
            "duration_s": t90 - t10 if t10 is not None and t90 is not None else None,
        }
    key = (180, .5)
    t10 = results["180_0.5"]["t10_clip_s"]
    t90 = results["180_0.5"]["t90_clip_s"]
    gain = .8 * nominal_final * PORE_VOLUME_UL
    delivered = PUMP_UL_S * (t90 - t10)
    report = {
        "source": "doi:10.18419/DARUS-4683", "license": "CC BY 4.0",
        "angle_deg": angle, "source_md5": spec["md5"], "frames": image.n_frames,
        "source_pixel_pitch_um": PITCH_UM, "frame_dt_s": dt,
        "wall_box_xyxy": [x0, y0, x1, y1], "interior_pixels": denominator,
        "assumed_uniform_depth_um": DEPTH_UM,
        "threshold_results": results,
        "reference_mask": {
            "contrast_threshold_counts": 180, "water_fraction_of_local_contrast": .5,
            "delta_saturation_by_definition": .8 * nominal_final,
            "pore_volume_gain_uL": gain, "nominal_pump_delivery_uL": delivered,
            "conditional_net_upstream_storage_release_uL": gain - delivered,
        },
        "not_established": [
            "complete thin-film detection", "independent 3D height profile",
            "camera-pressure hardware trigger", "curvature of 5.6 um fillet",
            "validated apparatus compliance", "capillary pressure prediction"
        ],
    }
    # Publish a small derived 10-Hz curve, not original images or TIFF.
    rows = [[angle, round(float(t[j]), 4),
             *[round(float(curves[k][j]), 5)
               for k in ((180, .4), (180, .5), (180, .6))]]
            for j in range(0, image.n_frames, 3)]
    return report, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--angle", type=int, required=True, choices=[60, 120])
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-tiff", type=Path)
    source.add_argument("--download-original", action="store_true")
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent / "outputs")
    args = parser.parse_args()
    raw = read_original(args.angle, args.input_tiff, args.download_original)
    report, rows = reconstruct(args.angle, raw)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / f"mask_{args.angle}_10hz.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["angle_deg", "clip_t_s", "mask_lo", "mask_reference", "mask_hi"])
        w.writerows(rows)
    (args.output_dir / f"mask_{args.angle}_audit.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print("No raw video saved or published; no camera-pressure alignment performed.")


if __name__ == "__main__":
    main()
