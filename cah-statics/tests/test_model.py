"""Cheap checks of the mathematics. These do not validate experimental accuracy."""
import numpy as np
from model import PlugModel, independent_minimization


def test():
    m = PlugModel(40)
    result = m.critical_and_post()
    before, after = result["before"], result["after"]
    np.testing.assert_allclose(result["critical_volume_uL"], .3216300097, rtol=2e-7)
    np.testing.assert_allclose(result["predicted_length_drop_um"], 209.6782165, rtol=2e-6)
    assert result["symmetric_fold_volume_uL"] < result["critical_volume_uL"]
    assert result["min_landing_eigenvalue"] > 0
    assert result["energy_drop"] > 0
    np.testing.assert_allclose(m.volume(*before[:2]), m.volume(*after[:2]), atol=1e-9)
    np.testing.assert_allclose(m.p(after[0],"left"), m.p(after[1],"right"), atol=1e-9)
    assert abs(m.constrained_stability(before)[0][0]) < 1e-8
    assert min(m.constrained_stability(m.symmetric_state(.76))[0]) > 0
    z = np.array([-.62,.73,.53,.75])
    analytic_e, analytic_a = m.gradients(z)
    for i in range(4):
        h = np.zeros(4)
        h[i] = 2e-6
        left, right = m.energy_area(z-h), m.energy_area(z+h)
        np.testing.assert_allclose(
            [(right[0]-left[0])/(4e-6), (right[1]-left[1])/(4e-6)],
            [analytic_e[i],analytic_a[i]], atol=3e-7)
    seed = m.symmetric_state(result["half_length_at_instability_mm"]-.0001)
    seed[:2] += .01
    seed[2] = m.theta(seed[0],"left")
    seed[3] = m.theta(seed[1],"right")
    landing = independent_minimization(m,result["critical_volume_uL"]-.0002,seed)
    np.testing.assert_allclose(landing[:2],after[:2],atol=.004)
    print("PASS mathematical consistency: gradients, curvature, full constrained stability and equal-volume landing.")
    print("NOT TESTED: experimental jump sizes, Young angle or measured wall profile.")


if __name__ == "__main__":
    test()
