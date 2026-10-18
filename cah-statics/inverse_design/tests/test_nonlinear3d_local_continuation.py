"""Offline validation of the geometry-fitted local nonlinear 3D solver.

Slow tests require CasADi's sparse Newton solver. No downloads, media,
GitHub Actions, or experimental data. Compared only with independent
neutral-cylinder analytic/Jacobi benchmarks, not observed bifurcations.
"""
from __future__ import annotations
import importlib.util
import math
import unittest
import numpy as np

_AVAILABLE=importlib.util.find_spec("casadi") is not None
if _AVAILABLE:
    from nonlinear3d_local_continuation import (
        physical_arclength_nodes,cylinder_point,local_continuation,
        constrained_hessian_min_eigenvalue,geometry_numpy)


@unittest.skipUnless(_AVAILABLE,"CasADi is optional")
class LocalThreeDTests(unittest.TestCase):
    states={-1:.7039537628915958,+1:.7088690566135175}
    exact_Q={-1:-.0873159482495,+1:-.0851924443800}
    exact_p={-1:-2.0387830031805305,+1:-2.0316750540586255}

    def test_physical_nodes_resolve_axial_distance(self):
        for phase,s in self.states.items():
            t=physical_arclength_nodes(129,phase,s,32)
            self.assertEqual(len(t),129)
            self.assertTrue(np.all(np.diff(t)>0))
            x=np.array([cylinder_point(float(z),phase,s,32)[0]
                        for z in t[64:]])
            dx=np.diff(x)
            self.assertLess((dx.max()-dx.min())/dx.mean(),1e-8)

    def test_exact_neutral_pressure_both_channels(self):
        for phase,s in self.states.items():
            r=local_continuation(129,3,phase,s,angles=(90.,))[0]
            self.assertLess(abs(r["pressure_signed_per_mm"]-self.exact_p[phase]),1e-4)
            self.assertLess(r["kkt_rms"],1e-8)
            self.assertLess(abs(r["volume_residual_mm3"]),1e-9)
            self.assertTrue(r["positive_orientation"])

    def test_nonlinear_susceptibility_approaches_exact(self):
        for phase,s in self.states.items():
            r=local_continuation(129,9,phase,s)
            self.assertLess(abs(r[1]["apparent_Q_mm2"]-self.exact_Q[phase]),.002)
            self.assertLess(r[1]["max_cylinder_displacement_um"],2)
            self.assertGreater(r[1]["min_transverse_jacobian"],.05)

    def test_small_epsilon_stability(self):
        for phase,s in self.states.items():
            r=local_continuation(65,5,phase,s,angles=(90.,89.95,89.9))
            self.assertLess(abs(r[1]["apparent_Q_mm2"]-r[2]["apparent_Q_mm2"]),2e-5)
            self.assertLess(abs(r[2]["delta_signed_pressure_per_mm"]
                                +2*math.cos(math.radians(89.9))/.38),5e-5)

    def test_constrained_hessian_positive_in_discrete_graph_class(self):
        for phase,s in self.states.items():
            r=constrained_hessian_min_eigenvalue(65,5,phase,s,angle_deg=89.9)
            self.assertGreater(r["min_projected_hessian_eigenvalue"],1e-4)

    def test_straight_channel_neutral_pressure(self):
        r=local_continuation(17,5,-1,.7,angles=(90.,89.9),straight=True)
        self.assertLess(abs(r[0]["pressure_signed_per_mm"]
                            +2*math.cos(math.radians(32))),1e-8)
        self.assertLess(r[0]["max_cylinder_displacement_um"],1e-6)
        self.assertLess(r[1]["apparent_Q_mm2"],0)

    def test_invalid_seed_and_even_mesh_rejected(self):
        with self.assertRaises(ValueError):
            physical_arclength_nodes(64,-1,self.states[-1],32)
        with self.assertRaises(ValueError):
            local_continuation(17,3,-1,self.states[-1],angles=(89.9,))


if __name__=="__main__":
    unittest.main()
