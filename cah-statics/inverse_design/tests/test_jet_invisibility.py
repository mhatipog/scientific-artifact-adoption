"""Regression tests of the fourth-jet non-identifiability construction."""
from __future__ import annotations
import math
import unittest
from spectral_isolation import Channel, critical_state
from jet_invisibility import (
    fourth_jet_sensitivity, invisible_flip, symmetric_bump)


class JetInvisibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g=Channel()
        cls.states=[critical_state(p,cls.g) for p in (-1,1)]

    def test_sensitivity_is_positive(self):
        for st in self.states:
            result=fourth_jet_sensitivity(st,self.g.theta_deg)
            self.assertGreater(result["H_prime_mm"],0)
            self.assertGreater(result["dB4_d_w4"],0)

    def test_exact_local_quartic_sign_reversal(self):
        for st in self.states:
            result=invisible_flip(st,self.g.theta_deg)
            self.assertAlmostEqual(result["modified_B4"],-st["B4"],places=9)
            self.assertLess(result["modified_B4"]*st["B4"],0)

    def test_critical_root_and_robin_data_unchanged_by_jet(self):
        # All input jets of p', kappa, L, and beta are of order <=3.
        for st in self.states:
            resp=invisible_flip(st,self.g.theta_deg)
            self.assertEqual(resp["critical_position_mm"],st["s_mm"])
            self.assertAlmostEqual(st["beta_wall"],st["beta_zero"],delta=2e-9)

    def test_bump_has_prescribed_fourth_derivative(self):
        for st in self.states:
            s=st["s_mm"]
            amplitude=invisible_flip(st,self.g.theta_deg)["delta_w4_per_mm3"]
            h=.004*1e-4
            value=symmetric_bump(s+h,s,.004,amplitude)
            self.assertAlmostEqual(24*value/h**4,amplitude,delta=abs(amplitude)*5e-8)
            self.assertEqual(symmetric_bump(s,s,.004,amplitude),0.)

    def test_even_under_reflection(self):
        st=self.states[0]
        s=st["s_mm"]
        for x in (-.712,-.709,-.703,-.699,0,.699,.703,.709,.712):
            self.assertAlmostEqual(symmetric_bump(x,s,.004,1e5),
                                   symmetric_bump(-x,s,.004,1e5),places=17)

    def test_local_width_perturbation_is_tiny_but_no_fabrication_claim(self):
        for st in self.states:
            amp=invisible_flip(st,self.g.theta_deg)["delta_w4_per_mm3"]
            s=st["s_mm"]
            values=[symmetric_bump(s+.004*j/1000,s,.004,amp)
                    for j in range(-1000,1001)]
            self.assertLess(max(map(abs,values)),2e-7) # below 0.2 nm in mm
            self.assertGreater(max(map(abs,values)),1e-8)

    def test_periodic_even_bump_and_overlap_guard(self):
        s=self.states[1]["s_mm"]
        for x in (.1,.703,.708,.712,1.04):
            y=symmetric_bump(x,s,.004,1e5,period_mm=.3)
            self.assertAlmostEqual(y,symmetric_bump(x+.3,s,.004,1e5,period_mm=.3),places=15)
            self.assertAlmostEqual(y,symmetric_bump(-x,s,.004,1e5,period_mm=.3),places=15)
        with self.assertRaises(ValueError):
            symmetric_bump(s,s,.004,1e5,period_mm=.006)


if __name__=="__main__":
    unittest.main()
