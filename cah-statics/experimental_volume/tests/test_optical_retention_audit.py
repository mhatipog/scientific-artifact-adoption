"""Offline regression tests for retained-pixel photometry. No network calls."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from optical_retention_audit import analyze, summarize, INTERIOR_RADII_PX, SETTINGS


def synthetic():
    first=np.full((512,512),1000,dtype=float)
    final=first.copy()
    before=first.copy()
    after=first.copy()
    final[90:380,125:430]=500
    # Stable-wet pixels progress from 0.64 to 0.94 optical units.
    before[150:200,150:200]=680
    after[150:200,150:200]=530
    # Stable dry control progresses only 0.05 units.
    before[250:290,250:290]=950
    after[250:290,250:290]=925
    # New-onset mask pixels progress 0.10 -> 0.80.
    before[210:240,210:240]=950
    after[210:240,210:240]=600
    return first,final,before,after,(116,69,438,391)


def row(rows,counts=180,fraction=.5,radius=7):
    hits=[r for r in rows if r["count_threshold"]==counts
          and r["fraction_threshold"]==fraction
          and r["retained_boundary_exclusion_px"]==radius]
    assert len(hits)==1
    return hits[0]


class PhotometricRetentionTests(unittest.TestCase):
    def test_summaries_and_empty(self):
        self.assertIsNone(summarize([])["delta_median"])
        self.assertEqual(summarize([-0.2,0,0.2])["count"],3)
        self.assertAlmostEqual(summarize([-0.2,0,0.2])["delta_median"],0)

    def test_mask_interiors_and_dry_control(self):
        r=analyze(*synthetic())
        self.assertEqual(len(r),len(SETTINGS)*len(INTERIOR_RADII_PX))
        m=row(r)
        self.assertGreater(m["retained_wet_count"],1000)
        self.assertAlmostEqual(m["retained_delta_median"],.3,delta=.005)
        self.assertAlmostEqual(m["stable_off_delta_median"],0,delta=1e-10)
        self.assertGreater(m["newly_wet_pixels"],800)
        self.assertLess(row(r,radius=15)["retained_wet_count"],
                        row(r,radius=3)["retained_wet_count"])

    def test_global_illumination_drift_is_removed(self):
        first,final,before,after,box=synthetic()
        expected=row(analyze(first,final,before,after,box))
        shifted=row(analyze(first,final,before+41.25,after-19.75,box))
        for key in ("retained_delta_median","stable_off_delta_median",
                    "newly_wet_pixels","retained_wet_count"):
            self.assertAlmostEqual(expected[key],shifted[key],delta=1e-10)

    def test_malformed_frame_fails(self):
        first,final,before,after,box=synthetic()
        with self.assertRaises(ValueError):
            analyze(first[:400],final,before,after,box)


if __name__=="__main__":
    unittest.main(verbosity=2)
