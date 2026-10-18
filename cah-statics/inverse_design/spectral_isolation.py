"""Neutral-plate continuum Jacobi-spectrum audit for Fourier-phase capillary pairs.

All lengths in mm. This is an exact spectrum of the local cylindrical
Young-Laplace Jacobi operator at a pressure-stationary meniscus, not a
non-neutral plate, full global branch, or experimental stability result.

Run: python spectral_isolation.py
Dependencies: scipy, numpy. No JAX, file writes, or GitHub Actions.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq


@dataclass(frozen=True)
class Channel:
    high_order: int = 5
    theta_deg: float = 32.0
    high_scale_mm: float = .10625
    width_mm: float = 1.0
    fundamental_mm: float = .05
    period_mm: float = .3
    depth_mm: float = .38

    def __post_init__(self):
        if self.high_order < 5 or self.high_order % 2 != 1:
            raise ValueError("High Fourier order must be odd and >=5")
        if not (0 < self.theta_deg < 90 and self.depth_mm > 0):
            raise ValueError("Out-of-range contact angle or depth")


def planar_jets(x: float, phase: int, g: Channel) -> dict:
    """Analytic p',p'',p''' and fixed-volume B4; no finite differences."""
    if phase not in (-1, 1):
        raise ValueError("Phase must be -1 or +1")
    k = 2 * math.pi / g.period_mm
    u = k * x
    m = g.high_order
    b = phase * g.high_scale_mm / m**4
    ds = []
    for n in range(5):
        dn = (g.width_mm if n == 0 else 0) + k**n * (
            g.fundamental_mm * math.cos(u + n * math.pi / 2)
            + b * (m**n * math.cos(m*u + n*math.pi/2)
                   - (m+2)**n * math.cos((m+2)*u+n*math.pi/2)))
        ds.append(dn)
    w, w1, w2, w3, w4 = ds
    z, z1, z2, z3 = w1/2, w2/2, w3/2, w4/2
    den = 1 + z*z
    theta_y = math.radians(g.theta_deg)
    theta = theta_y + math.atan(z)
    t1 = z1 / den
    t2 = z2/den - 2*z*z1*z1/den**2
    t3 = z3/den - (2*z1**3+6*z*z1*z2)/den**2 + 8*z**2*z1**3/den**3
    cs, sn = math.cos(theta), math.sin(theta)
    f0, f1, f2, f3 = (cs, -sn*t1, -cs*t1*t1-sn*t2,
                       sn*t1**3-3*cs*t1*t2-sn*t3)
    v0, v1, v2, v3 = (1/w, -w1/w**2,
                       2*w1*w1/w**3-w2/w**2,
                       -6*w1**3/w**4+6*w1*w2/w**3-w3/w**2)
    p1 = 2*(f1*v0 + f0*v1)
    p2 = 2*(f2*v0 + 2*f1*v1 + f0*v2)
    p3 = 2*(f3*v0 + 3*f2*v1 + 3*f1*v2 + f0*v3)
    coeff = (2*theta-math.pi)/8
    tan = math.tan(theta)
    sec2 = 1/cs**2
    h_f = coeff*sec2 + tan/4
    h_f1 = sec2/2 + 2*coeff*sec2*tan
    h_f2 = 1.5*sec2*tan + 4*coeff*sec2*tan*tan + 2*coeff*sec2**2
    h1 = w + 2*w*w1*h_f + w*w*h_f1*t1
    h2 = (w1 + 2*(w1*w1+w*w2)*h_f + 4*w*w1*h_f1*t1
          + w*w*(h_f2*t1*t1+h_f1*t2))
    kappa = 2*cs/w
    phi = math.pi/2-theta
    L = phi/kappa
    wall_curvature = w2/(2*den**1.5)
    # Free-contact-line Robin term including signed curvature of support.
    beta_wall = kappa/math.tan(theta_y) + wall_curvature/math.sin(theta_y)
    # At p'=0 the differentiated exact circular family is a Jacobi zero mode.
    beta_zero = kappa*math.tan(phi)
    return dict(w=w, w1=w1, w2=w2, p1=p1, p2=p2, p3=p3,
                B4=(h2*p2-h1*p3/3)/4, kappa=kappa, phi=phi, L=L,
                theta_app_deg=math.degrees(theta),
                wall_curvature=wall_curvature,
                beta_wall=beta_wall, beta_zero=beta_zero,
                jacobi_amplitude=1+w*w2/(4*den))


def critical_state(phase: int, g: Channel) -> dict:
    """Count local p' roots before selecting; never silently hop branches."""
    xs = np.linspace(.685, .715, 1001)
    vv = [planar_jets(float(x), phase, g)["p1"] for x in xs]
    ii = [j for j in range(len(xs)-1) if vv[j]*vv[j+1] < 0]
    if len(ii) != 1:
        raise ValueError(f"Expected one critical root, found {len(ii)}")
    s = brentq(lambda x: planar_jets(x, phase, g)["p1"],
               float(xs[ii[0]]), float(xs[ii[0]+1]), xtol=1e-14)
    d = planar_jets(s, phase, g)
    if (d["p2"] >= 0 or d["jacobi_amplitude"] <= 0
            or not (0 < d["phi"] < math.pi/2)):
        raise ValueError("Selected state fails the stated local stability gates")
    residual = d["beta_wall"]-d["beta_zero"]
    if abs(residual) > 2e-9:
        raise AssertionError("Physical wall-curvature Robin term disagrees with zero mode")
    return dict(phase=phase, s_mm=s, b_mm=phase*g.high_scale_mm/g.high_order**4,
                **d, robin_identity_residual_per_mm=residual)


def planar_eigenvalues(state: dict, count: int = 12) -> list[float]:
    """Exact Robin Sturm-Liouville eigenvalues, ordered with multiplicity."""
    L, kappa, beta = state["L"], state["kappa"], state["beta_wall"]
    bL = beta*L
    ev = []
    for j in range(count):
        even = brentq(lambda t: t*math.tan(t)-bL,
                      j*math.pi+1e-11, (j+.5)*math.pi-1e-11, xtol=1e-14)
        odd = brentq(lambda t: t/math.tan(t)+bL,
                     (j+.5)*math.pi+1e-11, (j+1)*math.pi-1e-11,
                     xtol=1e-14)
        ev.extend([(even/L)**2-kappa*kappa,
                   (odd/L)**2-kappa*kappa])
    return sorted(ev)[:count]


def full_3d_eigenvalues(state: dict, depth_mm: float, count: int = 8) -> list[float]:
    """At 90-degree plates, Neumann separation adds (n*pi/depth)^2."""
    ps = planar_eigenvalues(state, count + 2)
    values = sorted(lam+(n*math.pi/depth_mm)**2
                    for lam in ps for n in range(count+2))
    return values[:count]


def fe_planar_eigenvalues(state: dict, nodes: int) -> list[float]:
    """Independent conforming P1 FEM with physical Robin boundary term."""
    if nodes < 4:
        raise ValueError("At least four nodes")
    L, kappa, beta = state["L"], state["kappa"], state["beta_wall"]
    h = 2*L/(nodes-1)
    stiffness = np.zeros((nodes,nodes))
    mass = np.zeros((nodes,nodes))
    ek = np.array([[1., -1.],[-1., 1.]])/h
    em = h*np.array([[2., 1.],[1., 2.]])/6
    for i in range(nodes-1):
        stiffness[i:i+2, i:i+2] += ek
        mass[i:i+2, i:i+2] += em
    stiffness[0,0] += beta
    stiffness[-1,-1] += beta
    stiffness -= kappa*kappa*mass
    return [float(v) for v in eigh(stiffness, mass,
                                  subset_by_index=[0,1], driver="gvx")[0]]


def record(g: Channel) -> dict:
    minus, plus = [critical_state(phase, g) for phase in (-1, 1)]
    for s in (minus, plus):
        s["eigenvalues_first8_per_mm2"] = full_3d_eigenvalues(s,g.depth_mm)
        s["lateral_gap_per_mm2"] = planar_eigenvalues(s,2)[1]
        s["depth_gap_per_mm2"] = (math.pi/g.depth_mm)**2
        s["mode_crossover_depth_mm"] = math.pi/math.sqrt(s["lateral_gap_per_mm2"])
    gaps = (minus["lateral_gap_per_mm2"],plus["lateral_gap_per_mm2"])
    return dict(m=g.high_order, high_mode_magnitude_nm=1e6*g.high_scale_mm/g.high_order**4,
                minus=minus, plus=plus,
                lateral_gap_difference_per_mm2=abs(gaps[0]-gaps[1]),
                relative_gap_difference=2*abs(gaps[0]-gaps[1])/sum(gaps),
                max_first8_gap_difference_per_mm2=max(
                    abs(a-b) for a,b in zip(
                        minus["eigenvalues_first8_per_mm2"],
                        plus["eigenvalues_first8_per_mm2"])))


def report() -> dict:
    pair5 = record(Channel())
    asymptotic = [record(Channel(high_order=m, theta_deg=37.1451464361036,
                                 high_scale_mm=.046875))
                  for m in (15, 25, 51, 75, 101, 151, 201)]
    fem = {}
    for label,s in (("minus",pair5["minus"]),("plus",pair5["plus"])):
        exact = s["lateral_gap_per_mm2"]
        fem[label] = [
            dict(nodes=n, lambda0_per_mm2=vals[0],
                 lambda1_per_mm2=vals[1], gap_error_per_mm2=vals[1]-exact)
            for n in (33,65,129,257)
            for vals in (fe_planar_eigenvalues(s,n),)]
    return {
        "status":"local neutral-plate 3D continuum Jacobi spectrum; model only",
        "neutral_plate_angle_deg":90,
        "depth_mm":.38,
        "unit_eigenvalue":"mm^-2",
        "nominal_five_seven_pair":pair5,
        "fixed_degeneration_angle_high_order_sequence":asymptotic,
        "independent_robin_P1_FEM_convergence":fem,
        "limitations":[
            "Free-contact-line Young angle assumed, not independently measured",
            "No non-neutral plate 3D quartic or Hessian continuation",
            "Zero mode and gap apply only at the selected local critical cylinder",
            "No global energy minimality or earlier path instability proved",
            "Large-harmonic mathematical limits are not fabricable channel designs"]}


def csv_string(data: dict | None = None) -> str:
    """Stable, machine-readable rounded-free derived data for the SVG."""
    import csv
    import io
    if data is None:
        data = report()
    rows = ([data["nominal_five_seven_pair"]] +
            data["fixed_degeneration_angle_high_order_sequence"])
    fields = ("m", "high_mode_magnitude_nm", "B4_minus", "B4_plus",
              "s_minus_mm", "s_plus_mm",
              "gap_minus_per_mm2", "gap_plus_per_mm2",
              "lateral_gap_difference_per_mm2",
              "relative_gap_difference",
              "max_first8_gap_difference_per_mm2")
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for rec in rows:
        writer.writerow(dict(
            m=rec["m"], high_mode_magnitude_nm=rec["high_mode_magnitude_nm"],
            B4_minus=rec["minus"]["B4"], B4_plus=rec["plus"]["B4"],
            s_minus_mm=rec["minus"]["s_mm"], s_plus_mm=rec["plus"]["s_mm"],
            gap_minus_per_mm2=rec["minus"]["lateral_gap_per_mm2"],
            gap_plus_per_mm2=rec["plus"]["lateral_gap_per_mm2"],
            lateral_gap_difference_per_mm2=rec["lateral_gap_difference_per_mm2"],
            relative_gap_difference=rec["relative_gap_difference"],
            max_first8_gap_difference_per_mm2=rec["max_first8_gap_difference_per_mm2"]))
    return buffer.getvalue()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", action="store_true",
                        help="Print reproducible derived CSV to standard output")
    args = parser.parse_args()
    data = report()
    print(csv_string(data) if args.csv else json.dumps(data, indent=2), end="")

