"""Offline tests for exact neutral-plate plate-wetting response and FEM."""
from __future__ import annotations
import math
import unittest
from spectral_isolation import Channel, critical_state
from plate_susceptibility import (
    even_mode, susceptibility_exact, susceptibility_fem, p1_1d, csv_content, report)


class PlateSusceptibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g=Channel()
        cls.states={phase:critical_state(phase,cls.g) for phase in (-1,1)}

    def test_physical_wall_robin_identity(self):
        for s in self.states.values():
            self.assertLess(abs(s["beta_wall"]-s["kappa"]*math.tan(s["kappa"]*s["L"])),2e-9)
            self.assertGreater(s["beta_wall"],0)

    def test_even_mode_zero_and_excited(self):
        for s in self.states.values():
            k0,g0_sq=even_mode(0,s["L"],s["beta_wall"])
            k1,g1_sq=even_mode(1,s["L"],s["beta_wall"])
            self.assertAlmostEqual(k0,s["kappa"],delta=5e-9)
            self.assertGreater(k1,k0)
            self.assertGreater(g0_sq,0)
            self.assertGreater(g1_sq,0)

    def test_exact_response_phase_values(self):
        e0=susceptibility_exact(self.states[-1],.38,320)
        e1=susceptibility_exact(self.states[1],.38,320)
        self.assertAlmostEqual(e0["Q_mm2"],-.0873159482495,delta=3e-11)
        self.assertAlmostEqual(e1["Q_mm2"],-.0851924443800,delta=2e-11)
        self.assertGreater(e1["Q_mm2"]-e0["Q_mm2"],.0021)

    def test_spectral_tail_bound_against_longer_sum(self):
        for state in self.states.values():
            short=susceptibility_exact(state,.38,32)
            long=susceptibility_exact(state,.38,512)
            self.assertGreater(short["Q_mm2"],long["Q_mm2"])
            self.assertLess(short["Q_mm2"]-long["Q_mm2"],
                            short["spectral_tail_bound_mm2"])
            self.assertLess(long["spectral_tail_bound_mm2"],2e-11)

    def test_strict_negative_and_pressure_solvability(self):
        for state in self.states.values():
            r=susceptibility_exact(state,.38,48)
            self.assertLess(r["Q_mm2"],0)
            self.assertLess(r["higher_modes_mm2"],0)
            self.assertAlmostEqual(r["pressure_response_per_mm"],2/.38,places=12)

    def test_fem_homogeneous_pressure_constraint(self):
        for state in self.states.values():
            r=susceptibility_fem(state,.38,33,17)
            self.assertLess(r["stationarity_relative_residual"],1e-11)
            self.assertLess(abs(r["volume_residual_mm3"]),1e-11)
            self.assertAlmostEqual(r["multiplier_per_mm"],2/.38,delta=1e-8)

    def test_fem_second_order_convergence_both_phases(self):
        for state in self.states.values():
            exact=susceptibility_exact(state,.38,320)["Q_mm2"]
            errors=[]
            for nq,nz in ((17,9),(33,17),(65,33),(129,65)):
                f=susceptibility_fem(state,.38,nq,nz)
                errors.append(abs(f["Q_mm2"]-exact))
            self.assertLess(errors[-1],2.6e-5)
            self.assertTrue(all(a>b>0 for a,b in zip(errors[:-1],errors[1:])))
            self.assertTrue(all(3.8<a/b<4.2 for a,b in zip(errors[:-1],errors[1:])))

    def test_straight_channel_agrees_with_nonlinear_limit(self):
        angle=math.radians(32)
        kappa=2*math.cos(angle)
        L=(math.pi/2-angle)/kappa
        state=dict(L=L,kappa=kappa,beta_wall=kappa/math.tan(angle))
        exact=susceptibility_exact(state,.38,320)
        self.assertAlmostEqual(exact["Q_mm2"],-.0753987025,delta=5e-9)
        fe=susceptibility_fem(state,.38,129,65)
        self.assertLess(abs(fe["Q_mm2"]-exact["Q_mm2"]),6e-5)

    def test_depth_mesh_represents_actual_depth(self):
        _,M=p1_1d(9,.38/2)
        self.assertAlmostEqual(float(M.sum()),.38,places=12)

    def test_csv_is_deterministic_and_complete(self):
        x=report()
        csv=csv_content(x)
        self.assertEqual(len(csv.splitlines()),9)
        self.assertIn("absolute_error_mm2",csv.splitlines()[0])
        self.assertIn("-1,129,65",csv)
        self.assertIn("1,129,65",csv)

    def test_reject_invalid_inputs(self):
        state=self.states[-1]
        with self.assertRaises(ValueError):
            susceptibility_exact(state,0)
        with self.assertRaises(ValueError):
            susceptibility_exact(state,.38,1)
        bad={**state,"beta_wall":state["beta_wall"]+.01}
        with self.assertRaises(ValueError):
            susceptibility_exact(bad,.38,10)


if __name__=="__main__":
    unittest.main()
