"""Projection-only flow-budget audit of Yan et al. 95-degree source TIFF.

This module downloads no data on import. --download-ranges deliberately
requests only 206 byte ranges (header, first/final and seven burst frames)
from the original 2.54-GB TIFF. Never stores or republishes raw microscopy.
Credit: Yan et al. (2025), doi:10.18419/DARUS-4683, CC BY 4.0.

The optical class is NOT a calibrated depth-integrated liquid volume.
"""
from __future__ import annotations

import argparse
import csv
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from sparse_third_wettability import (
    BYTES_PER_FRAME, EXPECTED_DT, LAST, PITCH_UM, PUMP_UL_S,
    SOURCE_BYTES, SOURCE_ID, URL, find_pore, get_range, parse_first_ifd
)

FIRST_BURST = 4818
LAST_BURST = 4824
SETTINGS = ((100, 0.4), (180, 0.4), (180, 0.5), (180, 0.6), (260, 0.6))
WINDOWS = ((4820, 4822), (4822, 4824))
NOMINAL_MAX_DEPTH_UM = 120.0
TEST_INCREMENT_UM = 2.0


def accounting(added_px: int, lost_px: int, elapsed_s: float,
               pitch_um: float, flow_ul_s: float, max_depth_um: float,
               incremental_depth_um: float = TEST_INCREMENT_UM):
    """Conservation inequalities, with a maximal lost-pixel credit.

    Added pixels are initially BELOW the optical threshold. For the
    lower-bound redistribution result, incremental_depth_um must mean
    an independently verified LOWER BOUND on the AREA-AVERAGE increase
    in actual water-column thickness over newly positive pixels, not their final thickness.

    Lost positive pixels may release up to max_depth_um in total; other
    releases (unchanged-positive pixels, optical subthreshold liquid,
    apparatus upstream) count as unresolved redistribution/storage.
    """
    if not (added_px > lost_px >= 0 and elapsed_s > 0 and pitch_um > 0
            and flow_ul_s >= 0 and max_depth_um > 0 and
            incremental_depth_um >= 0):
        raise ValueError("Invalid mass-budget inputs / nonpositive net area")
    unit_area_mm2 = (pitch_um / 1000) ** 2
    added_mm2 = added_px * unit_area_mm2
    lost_mm2 = lost_px * unit_area_mm2
    pump_ul = flow_ul_s * elapsed_s  # 1 uL = 1 mm^3
    # Common, but assumes zero loss/thinning elsewhere in the mask.
    ideal_um = 1000 * pump_ul / (added_mm2 - lost_mm2)
    # Allow the disappearing pixels to drain their entire depth.
    loss_credit_ul = lost_mm2 * max_depth_um / 1000
    guarded_um = 1000 * (pump_ul + loss_credit_ul) / added_mm2
    # Unseen redistribution needed IF the prescribed incremental depth
    # is a certified lower bound on the newly positive area average (no device attribution).
    min_redistribution_ul = max(
        0.0, added_mm2 * incremental_depth_um / 1000
        - pump_ul - loss_credit_ul)
    return dict(added_area_mm2=added_mm2, lost_area_mm2=lost_mm2,
                pump_ul=pump_ul, ideal_critical_um=ideal_um,
                lost_pixel_guarded_critical_um=guarded_um,
                conditional_min_redistribution_ul=min_redistribution_ul)


def transitions(before: np.ndarray, after: np.ndarray):
    """Separate positive and negative transitions; net alone hides reversals."""
    before = np.asarray(before, dtype=bool)
    after = np.asarray(after, dtype=bool)
    if before.shape != after.shape:
        raise ValueError("Mask shape mismatch")
    return int(np.count_nonzero(after & ~before)), int(
        np.count_nonzero(before & ~after))


def classify(reference, final, current, roi, exterior, threshold, fraction):
    """Frozen Round-12 phase-contrast mask and exterior illumination correction."""
    reference = np.asarray(reference, float)
    final = np.asarray(final, float)
    current = np.asarray(current, float)
    roi = np.asarray(roi, bool)
    exterior = np.asarray(exterior, bool)
    if not (reference.shape == final.shape == current.shape ==
            roi.shape == exterior.shape) or not exterior.any():
        raise ValueError("Inconsistent images or no exterior pixels")
    contrast = reference - final
    drift = float(np.median((current-reference)[exterior]))
    return ((roi & (contrast > threshold)) &
            ((reference - current + drift) > fraction * contrast))


def analyze(first, final, burst, box):
    """Original arrays: first, last, and clip frames 4818..4824 inclusive."""
    if (first.shape != (512, 512) or final.shape != (512, 512) or
            burst.shape != (LAST_BURST-FIRST_BURST+1, 512, 512)):
        raise ValueError("Wrong original TIFF frame shape/count")
    x0,y0,x1,y1 = box
    roi = np.zeros((512,512),bool)
    roi[y0+3:y1-3, x0+3:x1-3] = True
    exterior = np.zeros_like(roi)
    exterior[:y0-10,180:330] = True
    exterior[y1+10:,180:330] = True
    ref = gaussian_filter(np.asarray(first,float), 1.)
    fin = gaussian_filter(np.asarray(final,float), 1.)
    blurred = {j:gaussian_filter(np.asarray(burst[j-FIRST_BURST],float),1.)
               for j in (4820,4822,4824)}
    rows=[]
    for i,j in WINDOWS:
        for counts,frac in SETTINGS:
            bef = classify(ref, fin, blurred[i], roi, exterior, counts, frac)
            aft = classify(ref, fin, blurred[j], roi, exterior, counts, frac)
            added,lost = transitions(bef,aft)
            values = accounting(added,lost,(j-i)*EXPECTED_DT,PITCH_UM,
                                PUMP_UL_S,NOMINAL_MAX_DEPTH_UM)
            rows.append(dict(frame_start=i,frame_end=j,
                count_threshold=counts, fraction_threshold=frac,
                added_pixels=added,lost_pixels=lost,
                net_pixels=added-lost,
                loss_to_add_ratio=lost/added, **values))
    return rows


def download_and_analyze():
    """Explicit opt-in to at most 4.8MB strictly bounded original ranges."""
    source=URL+str(SOURCE_ID)
    head,signed=get_range(source,0,4096,SOURCE_BYTES)
    tags,offset=parse_first_ifd(head)
    first,signed=get_range(signed,offset,BYTES_PER_FRAME,SOURCE_BYTES)
    final,signed=get_range(signed,offset+LAST*BYTES_PER_FRAME,
                           BYTES_PER_FRAME,SOURCE_BYTES)
    data,signed=get_range(signed,offset+FIRST_BURST*BYTES_PER_FRAME,
                          (LAST_BURST-FIRST_BURST+1)*BYTES_PER_FRAME,
                          SOURCE_BYTES)
    first=np.frombuffer(first,dtype=">u2").reshape(512,512)
    final=np.frombuffer(final,dtype=">u2").reshape(512,512)
    burst=np.frombuffer(data,dtype=">u2").reshape(-1,512,512)
    box=find_pore(first)
    if box != (116,69,438,391):
        raise ValueError("Source pore geometry differs from Round-12 audit")
    return analyze(first,final,burst,box)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-ranges",action="store_true",
        help="Opt in to ~4.8 MB of strict 206 reads; never a CI default")
    args=parser.parse_args()
    if not args.download_ranges:
        parser.error("Required explicit --download-ranges; no default download")
    data=download_and_analyze()
    writer=csv.DictWriter(sys.stdout,fieldnames=list(data[0]))
    writer.writeheader()
    writer.writerows(data)
    print("Original TIFF not saved; no camera-pressure alignment attempted.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
