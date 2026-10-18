"""Retained-pixel photometric audit of original Yan 95-degree bright-field TIFF.

WARNING: normalized brightness is not liquid volume or water-column depth.
Sparse original acquisition: DaRUS doi:10.18419/DARUS-4683, file 351440,
CC BY 4.0. No raw image data are committed or saved by this script.
Only the first/final frames and 4822/4824 are read via strict HTTP 206.
"""
from __future__ import annotations

import argparse
import csv
import sys

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter
from sparse_third_wettability import (
    BYTES_PER_FRAME, LAST, SOURCE_BYTES, SOURCE_ID, URL,
    find_pore, get_range, parse_first_ifd,
)

START_FRAME = 4822
END_FRAME = 4824
SETTINGS = ((100, .4), (180, .4), (180, .5), (180, .6), (260, .6))
INTERIOR_RADII_PX = (0, 3, 7, 15)


def summarize(values):
    """Describe a distribution without interpreting optical signal as volume."""
    a = np.asarray(values, dtype=float)
    if not len(a):
        return dict(count=0, delta_p10=None, delta_median=None, delta_p90=None,
                    fraction_delta_gt_0p1=None, fraction_delta_lt_neg_0p1=None)
    lo, mid, hi = np.percentile(a, (10, 50, 90))
    return dict(count=int(len(a)),
                delta_p10=float(lo), delta_median=float(mid),
                delta_p90=float(hi),
                fraction_delta_gt_0p1=float(np.mean(a > .1)),
                fraction_delta_lt_neg_0p1=float(np.mean(a < -.1)))


def analyze(first, final, before, after, box):
    """Frozen round-12 eligibility and masks, plus photometric phase response.

    The photographic numerator is reference-current+exterior_drift, and the
    denominator is the observed reference-final contrast at each pixel.
    This dimensionless ratio has no demonstrated relationship to column depth.
    The category masks remain frozen at the two selected frames. Eroded
    retained-wet masks reduce contamination by a moving binary boundary.
    """
    images=(first, final, before, after)
    if any(np.shape(a) != (512, 512) for a in images):
        raise ValueError("Expected exactly four 512x512 TIFF planes")
    x0,y0,x1,y1 = box
    if not (10 <= x0 < x1 <= 502 and 11 <= y0 < y1 <= 490):
        raise ValueError("Pore box outside original field of view")
    roi=np.zeros((512,512),bool)
    roi[y0+3:y1-3, x0+3:x1-3]=True
    exterior=np.zeros((512,512),bool)
    exterior[:y0-10,180:330]=True
    exterior[y1+10:,180:330]=True
    ref, fin, bef, aft=[gaussian_filter(np.asarray(a,float),1.) for a in images]
    contrast=ref-fin
    drift_before=float(np.median((bef-ref)[exterior]))
    drift_after=float(np.median((aft-ref)[exterior]))
    d_before=ref-bef+drift_before
    d_after=ref-aft+drift_after
    # Division only for positive, sufficiently contrasted eligible pixels.
    safe_denominator=np.where(contrast>0, contrast, 1.)
    progress_before=d_before/safe_denominator
    progress_after=d_after/safe_denominator
    delta=progress_after-progress_before
    rows=[]
    for counts,frac in SETTINGS:
        eligible=roi & (contrast>counts)
        wet_before=eligible & (d_before>frac*contrast)
        wet_after=eligible & (d_after>frac*contrast)
        retained=wet_before & wet_after
        added=~wet_before & wet_after
        removed=wet_before & ~wet_after
        stable_off=eligible & ~wet_before & ~wet_after
        if not np.any(retained) or not np.any(stable_off):
            raise ValueError("Original classification has no control group")
        control=summarize(delta[stable_off])
        distances=distance_transform_edt(retained)
        for radius in INTERIOR_RADII_PX:
            interior=retained & (distances>radius)
            stats=summarize(delta[interior])
            rows.append({
                "frame_start":START_FRAME, "frame_end":END_FRAME,
                "count_threshold":counts, "fraction_threshold":frac,
                "retained_boundary_exclusion_px":radius,
                "retained_wet_count":stats["count"],
                "retained_wet_before_median":
                    float(np.median(progress_before[interior]))
                    if interior.any() else None,
                "retained_wet_after_median":
                    float(np.median(progress_after[interior]))
                    if interior.any() else None,
                "retained_delta_p10":stats["delta_p10"],
                "retained_delta_median":stats["delta_median"],
                "retained_delta_p90":stats["delta_p90"],
                "retained_frac_delta_gt_0p1":stats["fraction_delta_gt_0p1"],
                "retained_frac_delta_lt_neg_0p1":
                    stats["fraction_delta_lt_neg_0p1"],
                "stable_off_count":control["count"],
                "stable_off_delta_median":control["delta_median"],
                "newly_wet_pixels":int(np.count_nonzero(added)),
                "reversed_wet_pixels":int(np.count_nonzero(removed)),
                "drift_before_gray":drift_before,
                "drift_after_gray":drift_after,
            })
    return rows


def download_and_analyze():
    url=URL+str(SOURCE_ID)
    header,signed=get_range(url,0,4096,SOURCE_BYTES)
    _tags,offset=parse_first_ifd(header)
    imgs=[]
    for frame in (0,LAST,START_FRAME,END_FRAME):
        data,signed=get_range(signed,offset+frame*BYTES_PER_FRAME,
                              BYTES_PER_FRAME,SOURCE_BYTES)
        imgs.append(np.frombuffer(data,dtype=">u2").reshape(512,512))
    box=find_pore(imgs[0])
    if box != (116,69,438,391):
        raise ValueError("Original pore boundary regression failed")
    return analyze(*imgs,box)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-ranges",action="store_true",
                        help="Manually permit four sparse source TIFF frame reads")
    args=parser.parse_args()
    if not args.download_ranges:
        parser.error("Explicit --download-ranges required, no CI default")
    rows=download_and_analyze()
    writer=csv.DictWriter(sys.stdout,fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    print("Derived optical responses only. No physical volume inferred.",
          file=sys.stderr)


if __name__=="__main__":
    main()
