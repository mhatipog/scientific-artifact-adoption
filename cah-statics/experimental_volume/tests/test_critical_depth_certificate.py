"""Network-free regression tests for the conditional areal flow budget."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from critical_depth_certificate import accounting, analyze, classify, transitions


class ProjectionBudgetTests(unittest.TestCase):
    def test_zero_reversal_reduces_to_pump_over_new_area(self):
        r = accounting(4, 0, 1, 1000, 1, 1000, 500)
        self.assertAlmostEqual(r["ideal_critical_um"], 250)
        self.assertAlmostEqual(r["lost_pixel_guarded_critical_um"], 250)
        self.assertAlmostEqual(r["conditional_min_redistribution_ul"], 1)

    def test_maximal_full_depth_reversal_credit(self):
        r = accounting(3, 1, 1, 1000, 1, 1000, 1000)
        self.assertAlmostEqual(r["ideal_critical_um"], 500)
        self.assertAlmostEqual(r["lost_pixel_guarded_critical_um"], 2000/3)
        self.assertAlmostEqual(r["conditional_min_redistribution_ul"], 1)

    def test_activation_and_reversal_are_not_net_count(self):
        before = np.array([[0, 1, 0], [1, 0, 0]], dtype=bool)
        after = np.array([[1, 1, 1], [0, 0, 0]], dtype=bool)
        self.assertEqual(transitions(before, after), (2, 1))
        with self.assertRaises(ValueError):
            transitions(before, after[:, 1:])

    def test_exterior_illumination_correction(self):
        ref=np.full((2,3),100.)
        final=np.zeros((2,3))
        cur=np.array([[50,50,110],[100,100,110]],float)
        roi=np.array([[1,1,0],[0,0,0]],bool)
        exterior=np.array([[0,0,1],[0,0,1]],bool)
        mask=classify(ref,final,cur,roi,exterior,1,.55)
        self.assertEqual(mask.tolist(),[[True,True,False],[False,False,False]])

    def test_guards(self):
        with self.assertRaises(ValueError):
            accounting(1, 1, 1, 1, 1, 100)
        with self.assertRaises(ValueError):
            accounting(1, 0, 0, 1, 1, 100)

    def test_synthetic_increasing_front_no_network(self):
        first=np.full((512,512),2000.,float)
        final=first.copy()
        final[72:388,119:435]=1000
        frames=np.repeat(first[None,:,:],7,axis=0)
        for idx,width in ((2,10),(4,20),(6,40)):
            frames[idx,72:388,119:119+width]=1000
        rows=analyze(first,final,frames,(116,69,438,391))
        self.assertEqual(len(rows),10)
        self.assertTrue(all(r["net_pixels"]>0 for r in rows))

    def test_conditional_depth_not_claimed_as_measured(self):
        a=accounting(100,1,.1,10,.00001,120,2)
        # Greater measured incremental thickness would require more
        # non-pump redistribution; the calculation is conditional.
        b=accounting(100,1,.1,10,.00001,120,5)
        self.assertGreater(b["conditional_min_redistribution_ul"],
                           a["conditional_min_redistribution_ul"])


if __name__ == "__main__":
    unittest.main()
