"""Offline algebraic and falsifiable derivative-sign tests (no network)."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from phase_pair_window import Geometry, jets, stationary_roots
from pressure_volume_bridge import (
    bridge_at_stationary, pressure_third_derivative_sign_bound,
)


class PressureVolumeBridgeTests(unittest.TestCase):
    def test_general_chain_rule_polynomial_not_specific_geometry(self):
        # H'(s)=A,H''(s)=B at a stationary single-front p(Q).
        for A,B,pQQ,pQQQ in ((.2,-3,-2,15),(3,24,7,-19),(6,-44,-5,0)):
            p1=0
            p2=pQQ*A*A
            p3=pQQQ*A**3+3*pQQ*A*B
            r=bridge_at_stationary(A,B,p1,p2,p3)
            self.assertAlmostEqual(r.pressure_curvature_in_Q,pQQ,places=11)
            self.assertAlmostEqual(r.pressure_third_derivative_in_Q,pQQQ,
                                   places=10)
            self.assertAlmostEqual(r.identity_residual,0,delta=1e-9)

    def test_both_verified_phase_geometries(self):
        for sign,expected_p3,expected_b4 in (
                (-1,109.22174126199995,-2413.486980993566),
                (+1,-534.8302697928979,3605.9306784584574)):
            g=Geometry(b=sign*.00017)
            x=stationary_roots(g)[0]
            d=jets(x,g)
            r=bridge_at_stationary(
                d["H1"],d["H2"],d["p1"],d["p2"],d["p3"])
            self.assertAlmostEqual(r.pressure_third_derivative_in_Q,
                                   expected_p3,delta=1e-6)
            self.assertAlmostEqual(r.quartic_from_geometry,expected_b4,
                                   delta=1e-6)
            self.assertAlmostEqual(r.quartic_from_volume,expected_b4,
                                   delta=1e-9)
            self.assertLess(abs(r.identity_residual),1e-9)

    def test_interval_bound_on_known_cubic_quintic(self):
        # Synthetic analytic function, NOT experimental pressures.
        gamma=.001
        h=.02
        f=lambda q:q**3+gamma*q**5
        cert=pressure_third_derivative_sign_bound(
            f(-2*h),f(-h),f(h),f(2*h),h,1e-7,120*gamma)
        self.assertEqual(cert.third_derivative_sign,"positive")
        self.assertEqual(cert.local_quartic_sign,"negative")
        self.assertLess(abs(cert.third_derivative_estimate-6),
                        cert.worst_case_error_bound)
        # A nominal negative-p''' curve reverses the local quartic sign.
        negative=pressure_third_derivative_sign_bound(
            -f(-2*h),-f(-h),-f(h),-f(2*h),h,1e-7,120*gamma)
        self.assertEqual(negative.local_quartic_sign,"positive")

    def test_no_unjustified_certificate_with_large_noise(self):
        h=.02
        cert=pressure_third_derivative_sign_bound(
            -8*h**3,-h**3,h**3,8*h**3,h,.1,0)
        self.assertEqual(cert.third_derivative_sign,"indeterminate")

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            bridge_at_stationary(0,0,0,1,1)
        with self.assertRaises(ValueError):
            bridge_at_stationary(1,0,.1,1,1)
        with self.assertRaises(ValueError):
            pressure_third_derivative_sign_bound(1,2,3,4,0,.01,5)
        with self.assertRaises(ValueError):
            pressure_third_derivative_sign_bound(1,2,3,4,1,-.01,5)


if __name__=="__main__":
    unittest.main(verbosity=2)
