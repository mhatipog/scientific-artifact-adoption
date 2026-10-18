"""Offline regression, Robin boundary audit and P1 mesh tests, no JAX."""
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spectral_isolation import (
    Channel, critical_state, planar_jets, planar_eigenvalues,
    full_3d_eigenvalues, fe_planar_eigenvalues, record
)
from phase_pair_window import Geometry, jets, stationary_roots


class SpectralIsolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = record(Channel())

    def test_independent_analytic_jets_match_existing_model(self):
        for phase in (-1,1):
            s = self.base["minus" if phase < 0 else "plus"]["s_mm"]
            old = jets(s,Geometry(b=phase*.00017))
            new = planar_jets(s,phase,Channel())
            for k in ("p1","p2","p3","B4","w","w1","w2"):
                self.assertAlmostEqual(old[k],new[k],delta=1e-6)

    def test_existing_critical_roots_and_quartics(self):
        self.assertAlmostEqual(self.base["minus"]["s_mm"],.7039537628915958,delta=1e-9)
        self.assertAlmostEqual(self.base["plus"]["s_mm"],.7088690566135175,delta=1e-9)
        self.assertAlmostEqual(self.base["minus"]["B4"],-2413.48698099,delta=.001)
        self.assertAlmostEqual(self.base["plus"]["B4"],3605.93067846,delta=.001)

    def test_physical_robin_includes_support_curvature(self):
        for s in (self.base["minus"],self.base["plus"]):
            self.assertLess(abs(s["robin_identity_residual_per_mm"]),2e-9)
            theta = math.radians(32)
            beta_flat = s["kappa"]/math.tan(theta)
            self.assertGreater(abs(s["beta_wall"]-beta_flat),5)
            self.assertGreater(s["jacobi_amplitude"],0)

    def test_robin_zero_mode_only_at_pressure_stationarity(self):
        for phase in (-1,1):
            root = self.base["minus" if phase < 0 else "plus"]["s_mm"]
            for dx in (-.001,.001):
                d = planar_jets(root+dx,phase,Channel())
                self.assertGreater(abs(d["p1"]),.1)
                self.assertGreater(abs(d["beta_wall"]-d["beta_zero"]),.4)

    def test_exact_continuum_neutral_spectrum_and_depth_separation(self):
        for phase in ("minus","plus"):
            s = self.base[phase]
            planar = planar_eigenvalues(s,6)
            full = full_3d_eigenvalues(s,.38,8)
            self.assertLess(abs(planar[0]),1e-9)
            self.assertTrue(all(v > 0 for v in planar[1:]))
            self.assertTrue(all(full[i] <= full[i+1] for i in range(7)))
            self.assertAlmostEqual(full[1],planar[1],delta=1e-9)
            self.assertAlmostEqual(full[4],(math.pi/.38)**2,delta=1e-8)
            self.assertGreater(s["mode_crossover_depth_mm"],.88)
            self.assertLess(s["mode_crossover_depth_mm"],.90)

    def test_conforming_FEM_second_order_convergence(self):
        for phase in ("minus","plus"):
            s = self.base[phase]
            exact = s["lateral_gap_per_mm2"]
            errs = []
            nulls = []
            for n in (33,65,129,257):
                v0,v1 = fe_planar_eigenvalues(s,n)
                self.assertGreater(v0,0)
                self.assertGreater(v1,exact)
                nulls.append(v0)
                errs.append(v1-exact)
            for i in range(3):
                self.assertGreater(errs[i]/errs[i+1],3.8)
                self.assertLess(errs[i]/errs[i+1],4.2)
                self.assertGreater(nulls[i]/nulls[i+1],3.8)
                self.assertLess(nulls[i]/nulls[i+1],4.2)
            self.assertLess(errs[-1],.00071)
            self.assertLess(nulls[-1],.00005)

    def test_nominal_mismatch_quartic_reversal_but_small_gap_change(self):
        g0,g1 = (self.base[phase]["lateral_gap_per_mm2"]
                 for phase in ("minus","plus"))
        self.assertLess(self.base["minus"]["B4"],0)
        self.assertGreater(self.base["plus"]["B4"],0)
        self.assertLess(2*abs(g0-g1)/(g0+g1),.0015)
        self.assertGreater(min(g0,g1),12)

    def test_asymptotic_reversal_with_closing_eigenvalue_gap(self):
        by_order = [record(Channel(high_order=m,theta_deg=37.1451464361036,
                                   high_scale_mm=.046875))
                    for m in (15,51,201)]
        for rec in by_order:
            self.assertLess(rec["minus"]["B4"]*rec["plus"]["B4"],0)
            self.assertLess(rec["lateral_gap_difference_per_mm2"]*
                            rec["m"]**2,2.)
            self.assertLess(abs(rec["minus"]["robin_identity_residual_per_mm"]),2e-9)
            self.assertLess(abs(rec["plus"]["robin_identity_residual_per_mm"]),2e-9)
        diffs = [rec["max_first8_gap_difference_per_mm2"] for rec in by_order]
        self.assertGreater(diffs[0],diffs[1])
        self.assertGreater(diffs[1],diffs[2])
        self.assertLess(diffs[-1],.0005)
        self.assertGreater(min(abs(by_order[-1][p]["B4"]) for p in ("minus","plus")),400)

    def test_csv_and_svg_artifact_reproducers(self):
        import csv
        import io
        import xml.etree.ElementTree as ET
        from plot_spectral_isolation import svg_string
        from spectral_isolation import csv_string
        csv_data = csv_string()
        rows = list(csv.DictReader(io.StringIO(csv_data)))
        self.assertEqual([int(r["m"]) for r in rows],
                         [5,15,25,51,75,101,151,201])
        svg = svg_string()
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))
        self.assertEqual(svg.count("<polyline"),3)
        root_folder = Path(__file__).resolve().parents[1]
        for path,expected in (
                (root_folder/"derived"/"neutral_plate_jacobi_spectrum.csv",csv_data),
                (root_folder/"figures"/"jacobi_spectral_isolation.svg",svg)):
            if path.exists():
                self.assertEqual(path.read_text(),expected)

    def test_invalid_input_rejected(self):
        with self.assertRaises(ValueError): Channel(high_order=4)
        with self.assertRaises(ValueError): Channel(theta_deg=100)
        with self.assertRaises(ValueError): planar_jets(.7,0,Channel())
        with self.assertRaises(ValueError): fe_planar_eigenvalues(self.base["minus"],2)


if __name__=="__main__":
    unittest.main()
