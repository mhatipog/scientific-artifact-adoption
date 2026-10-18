"""Offline regression of full reduced-model pressure equalities to the normal form."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from local_branch_validation import (
    DISPLACEMENTS_MM, PHASE_SIGNS, branch_point, rows,
)


class LocalBranchValidationTests(unittest.TestCase):
    def test_exact_pressure_balance_and_sign_of_branch_side(self):
        data=rows()
        self.assertEqual(len(data),8)
        for r in data:
            self.assertLess(abs(r["pressure_equality_residual"]),2e-12)
            self.assertEqual(r["mean_volume_offset_mm2"]>0,
                             r["phase_sign"]==-1)
            self.assertGreater(r["half_volume_imbalance_mm2"],0)

    def test_second_order_convergence_of_volume_normal_form(self):
        # Regression tolerances are for this chosen reduced numerical model,
        # not preregistered experimental acceptance thresholds.
        for sign in PHASE_SIGNS:
            err=[branch_point(sign,c)["relative_slope_error"]
                 for c in DISPLACEMENTS_MM]
            self.assertTrue(all(err[i]>err[i+1] for i in range(3)))
            self.assertLess(err[0],.05)
            self.assertLess(err[-1],.001)
            # Halving c should asymptotically reduce O(c²) error ~fourfold.
            self.assertGreater(err[2]/err[3],3.5)
            self.assertLess(err[2]/err[3],4.5)

    def test_invalid_inputs(self):
        for sign,c in ((0,.001),(1,0),(1,.02)):
            with self.assertRaises(ValueError):
                branch_point(sign,c)


if __name__=="__main__":
    unittest.main(verbosity=2)
