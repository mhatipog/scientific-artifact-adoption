"""Offline exact-mode, derivative, fold and local-geometry checks."""
import math
import unittest
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dataclasses import replace

from phase_pair_window import (
    Geometry, clearance, exact_low_order_correlations, jets, pair,
    stationary_fold_creation, stationary_roots, stress_grid,
    width_separation_um
)


class PhasePairWindowTests(unittest.TestCase):
    def test_integer_fourier_correlation_parity(self):
        counts=exact_low_order_correlations(5)
        self.assertEqual(len(counts),5)
        self.assertTrue(all(c>0 for c in counts))
        with self.assertRaises(ValueError):
            exact_low_order_correlations(6)

    def test_previously_reported_75nm_numerical_baseline(self):
        d=pair(Geometry(theta_deg=36),.000075,.000075)
        self.assertTrue(d["ok"])
        self.assertAlmostEqual(d["B4_minus"],-1505.931,delta=.2)
        self.assertAlmostEqual(d["B4_plus"],812.389,delta=.2)
        self.assertAlmostEqual(d["s_minus_mm"],.70019156,delta=1e-6)
        self.assertAlmostEqual(d["s_plus_mm"],.69945579,delta=1e-6)

    def test_analytic_jets_against_independent_differencing(self):
        g=Geometry(theta_deg=32,b=.00017)
        x=stationary_roots(g)[0]
        analytic=jets(x,g)
        for step in (1e-6,1e-7):
            lo,hi=jets(x-step,g),jets(x+step,g)
            for original,antiderivative in (
                    ("p2","p1"),("p3","p2"),("H2","H1")):
                numerical=(hi[antiderivative]-lo[antiderivative])/(2*step)
                self.assertLess(abs(numerical-analytic[original])/
                                max(1,abs(analytic[original])),2e-6)

    def test_exact_nominal_extrema_and_width_difference(self):
        # b<a/24 preserves base sinusoid's strict global extremes.
        g=Geometry()
        self.assertLess(g.b,g.a/24)
        for sign in (-1,1):
            geo=replace(g,b=sign*g.b)
            k=2*math.pi/geo.period
            for x,expect in ((0,1.05),(math.pi/k,.95)):
                width=geo.nominal_width+geo.a*math.cos(k*x)+geo.b*(
                    math.cos(5*k*x)-math.cos(7*k*x))
                self.assertAlmostEqual(width,expect,places=12)
        self.assertAlmostEqual(width_separation_um(g.b),.6574657,delta=1e-6)

    def test_nominal_local_onsets_and_neutral_jacobi(self):
        r=pair(Geometry(),.00017,.00017)
        self.assertTrue(r["ok"])
        self.assertLess(r["B4_minus"],-2000)
        self.assertGreater(r["B4_plus"],3000)
        self.assertGreater(r["jacobi_minus"],3)
        self.assertGreater(r["jacobi_plus"],3)
        self.assertGreater(r["clearance_minus_mm"],0)
        self.assertGreater(r["clearance_plus_mm"],0)

    def test_405_discrete_uncertainty_cases(self):
        rows=stress_grid()
        self.assertEqual(len(rows),405)
        self.assertTrue(all(r["ok"] and
            r["clearance_minus_mm"]>0 and r["clearance_plus_mm"]>0
            for r in rows))
        self.assertGreater(min(min(-r["B4_minus"],r["B4_plus"])
                               for r in rows),200)
        self.assertGreater(min(min(-r["p2_minus"],-r["p2_plus"])
                               for r in rows),40)
        self.assertGreater(min(min(r["jacobi_minus"],r["jacobi_plus"])
                               for r in rows),2.7)

    def test_competing_pressure_fold_born_at_finite_amplitude(self):
        p=stationary_fold_creation()
        self.assertAlmostEqual(p["b_mm"]*1000,.241991386,delta=5e-7)
        self.assertLess(p["residual_norm"],1e-7)
        g=Geometry(a=.049,period=.303,theta_deg=34)
        before=stationary_roots(replace(g,b=p["b_mm"]-.000002))
        after=stationary_roots(replace(g,b=p["b_mm"]+.000002))
        self.assertEqual(len(before),1)
        self.assertEqual(len(after),3)


if __name__=="__main__":
    unittest.main()
