"""Exact local transfer from single-front p(q) to an additive two-front plug.

For q = H(x), single-front pressure p(x), additive constrained energy
E(q1)+E(q2), E'(q) = -p(q), q1+q2 = 2q0:
  quadratic in volume imbalance delta: -p_q'(q0) delta**2
  quartic at p_q'(q0)=0: -p_q'''(q0) delta**4/12.
At x=s and center shift c, delta=H'(s)c+O(c**3),
  B4 = [3 H'' p'' - H' p''']/12 = -(H')**4 p_q'''/12.
These are IDENTITIES conditional on additivity and other modes' stability,
not global stability results or an original experimental discovery.
"""
from __future__ import annotations

from math import isclose


def quartic_from_geometry(h1, h2, p2, p3):
    if h1 <= 0:
        raise ValueError("Volume coordinate must increase with contact-line position")
    return (3.0 * h2 * p2 - h1 * p3) / 12.0


def single_front_third_q_derivative_at_stationary(h1, h2, p2, p3):
    if h1 <= 0:
        raise ValueError("Non-invertible volume coordinate")
    # Chain rule at p_x'=0; all derivatives on RHS are with respect to x.
    return p3 / h1**3 - 3 * p2 * h2 / h1**4


def transfer_coefficient(h1, h2, p2, p3):
    pq3 = single_front_third_q_derivative_at_stationary(h1, h2, p2, p3)
    return -(h1**4) * pq3 / 12


def synthetic_energy(c, a, b, alpha, beta):
    """Exactly integrate a polynomial model; solve exact-volume constraint."""
    def h(x):
        return x + a*x*x/2 + b*x*x*x/6

    def hd(x):
        return 1 + a*x + b*x*x/2

    # p(x)=1+alpha*x²/2+beta*x³/6; E(x)=-integral[p(x)H'(x)] dx.
    p = [1., 0., alpha/2, beta/6]
    hp = [1., a, b/2]
    product = [0.] * (len(p)+len(hp)-1)
    for i, pi in enumerate(p):
        for j, hj in enumerate(hp):
            product[i+j] += pi*hj

    def energy(x):
        return -sum(coef*x**(i+1)/(i+1) for i,coef in enumerate(product))

    s = 0.
    for _ in range(12):
        s -= (h(s+c)+h(s-c))/(hd(s+c)+hd(s-c))
    assert abs(h(s+c)+h(s-c)) < 1e-12
    return energy(s+c)+energy(s-c)


def run_checks():
    # Exact algebraic identity and convergent fixed-Q finite perturbation test.
    for beta in (-0.4, 1.8):
        b4 = quartic_from_geometry(1., .2, 2., beta)
        assert isclose(b4, transfer_coefficient(1., .2, 2., beta),
                       abs_tol=1e-14)
        near = synthetic_energy(.025, .2, .1, 2., beta) / .025**4
        assert abs(near-b4) < 4e-6, (beta, near, b4)
        assert (b4 > 0) == (beta < 1.2)
    # At a degenerate p_q'=p_q'''=0 the first possible even term is sixth.
    # E'(q)=-p0 - a*q^5/120 gives E(q)=-p0*q-a*q^6/720.
    # E(delta)+E(-delta)-2E(0)=-a*delta^6/360 exactly.
    a,delta = 2.5, .2
    e = lambda q: -1.2*q-a*q**6/720
    assert isclose(e(delta)+e(-delta),-a*delta**6/360,abs_tol=1e-16)
    print("PASS exact pressure-volume transfer, fixed-Q quartic, both signs, sixth-order degeneracy")


if __name__ == "__main__":
    run_checks()
