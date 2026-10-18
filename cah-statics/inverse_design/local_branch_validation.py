"""Direct nonlinear pressure-balance branch check for the two-front reduction.

For fixed center displacement c, solve p(s+c)=p(s-c) for s near the
pressure stationary point. Integrate H'(x) to obtain the *actual*
constrained mean-front volume q and imbalance eta on this local branch.
Compare eta^2/(q-qc) with -6 P''(qc)/P'''(qc). No 3D or experimental claim.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys

from scipy.integrate import quad
from scipy.optimize import brentq
from phase_pair_window import Geometry, jets, stationary_roots
from pressure_volume_bridge import bridge_at_stationary

DISPLACEMENTS_MM=(.004,.002,.001,.0005)
PHASE_SIGNS=(-1,1)


def pressure(x,g):
    d=jets(x,g)
    return 2*math.cos(d["apparent_angle_rad"])/d["w"]


def branch_point(sign,c):
    if sign not in PHASE_SIGNS or not (0<c<=.004):
        raise ValueError("Only declared phase pair and local displacements")
    g=Geometry(b=sign*.00017)
    roots=stationary_roots(g)
    if len(roots)!=1:
        raise ValueError("Local pressure stationary root not isolated")
    xc=roots[0]
    d=jets(xc,g)
    v=bridge_at_stationary(
        d["H1"],d["H2"],d["p1"],d["p2"],d["p3"])
    normal_slope=(-6*v.pressure_curvature_in_Q/
                  v.pressure_third_derivative_in_Q)
    def balance(s):
        return pressure(s+c,g)-pressure(s-c,g)
    s=brentq(balance,xc-.025,xc+.025,xtol=1e-14)
    # The reference integration origin is xc. Arbitrary H integration
    # constants cancel from both q-qc and the front-volume imbalance.
    q1=quad(lambda x:jets(x,g)["H1"],xc,s+c,epsabs=1e-12)[0]
    q2=quad(lambda x:jets(x,g)["H1"],xc,s-c,epsabs=1e-12)[0]
    q_shift=(q1+q2)/2
    eta=(q1-q2)/2
    if abs(q_shift)<1e-12:
        raise ValueError("Branch too near degeneracy for a stable ratio")
    numerical_slope=eta*eta/q_shift
    return dict(phase_sign=sign,displacement_mm=c,
                stationary_x_mm=xc,nonlinear_s_mm=s,
                mean_volume_offset_mm2=q_shift,
                half_volume_imbalance_mm2=eta,
                derivative_normal_form_slope_mm2=normal_slope,
                full_balance_slope_mm2=numerical_slope,
                relative_slope_error=abs(numerical_slope/normal_slope-1),
                pressure_equality_residual=balance(s))


def rows():
    return [branch_point(sign,c) for sign in PHASE_SIGNS
            for c in DISPLACEMENTS_MM]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv",action="store_true")
    args=parser.parse_args()
    data=rows()
    if args.csv:
        writer=csv.DictWriter(sys.stdout,fieldnames=list(data[0]))
        writer.writeheader();writer.writerows(data)
    else:
        for r in data:
            print(r)


if __name__=="__main__":
    main()
