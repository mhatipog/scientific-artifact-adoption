"""Numerical audit of a more resolvable 5/7 Fourier-phase capillary pair.

A *planar*, local, fixed-volume two-meniscus normal-form test with an
exact neutral-horizontal-plate 3D cylindrical stationary lift. NOT a
validated non-neutral 3D calculation, full-Hessian continuation, or
fabricated experiment. Reimplements capillary derivatives analytically,
independent of mpmath finite differencing used in Round 7/8.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import csv
import json
import math
import sys

import numpy as np
from scipy.optimize import brentq, root as nonlinear_root


@dataclass(frozen=True)
class Geometry:
    a: float = .05             # fundamental width amplitude, mm
    period: float = .3         # mm
    theta_deg: float = 32.     # UNMEASURED Young sidewall angle
    b: float = .00017          # signed higher-mode width amplitude, mm
    nominal_width: float = 1.  # mm


def jets(x: float, g: Geometry):
    """Exact derivatives p',p'',p''',H',H'' of Round-2 planar formula."""
    k = 2*math.pi/g.period
    u = k*x
    w = g.nominal_width + g.a*math.cos(u) + g.b*(
        math.cos(5*u)-math.cos(7*u))
    derivatives = [w]
    for order in range(1, 5):
        d = k**order*(g.a*math.cos(u+order*math.pi/2) + g.b*(
            5**order*math.cos(5*u+order*math.pi/2) -
            7**order*math.cos(7*u+order*math.pi/2)))
        derivatives.append(d)
    w,w1,w2,w3,w4 = derivatives
    z,z1,z2,z3 = w1/2,w2/2,w3/2,w4/2
    den = 1+z*z
    theta = math.radians(g.theta_deg)+math.atan(z)
    t1 = z1/den
    t2 = z2/den-2*z*z1*z1/den**2
    t3 = z3/den-(2*z1**3+6*z*z1*z2)/den**2+8*z*z*z1**3/den**3
    cs, sn = math.cos(theta),math.sin(theta)
    gv,gv1 = cs,-sn*t1
    gv2 = -cs*t1*t1-sn*t2
    gv3 = sn*t1**3-3*cs*t1*t2-sn*t3
    inv=1/w
    inv1=-w1/w**2
    inv2=2*w1*w1/w**3-w2/w**2
    inv3=-6*w1**3/w**4+6*w1*w2/w**3-w3/w**2
    p1=2*(gv1*inv+gv*inv1)
    p2=2*(gv2*inv+2*gv1*inv1+gv*inv2)
    p3=2*(gv3*inv+3*gv2*inv1+3*gv1*inv2+gv*inv3)
    sec2=1/cs**2
    tan=math.tan(theta)
    coeff=(2*theta-math.pi)/8
    f=coeff*sec2+tan/4
    f1=sec2/2+2*coeff*sec2*tan
    f2=1.5*sec2*tan+4*coeff*sec2*tan*tan+2*coeff*sec2**2
    H1=w+2*w*w1*f+w*w*f1*t1
    H2=w1+2*(w1*w1+w*w2)*f+4*w*w1*f1*t1 + w*w*(
        f2*t1*t1+f1*t2)
    B4=(H2*p2-H1*p3/3)/4
    jacobi=1+w*w2/(4*(1+(w1/2)**2))
    return dict(p1=p1,p2=p2,p3=p3,H1=H1,H2=H2,B4=B4,
                w=w,w1=w1,w2=w2,apparent_angle_rad=theta,
                jacobi_amplitude=jacobi)


def stationary_roots(g: Geometry, samples: int = 451):
    """Isolate sign-change p' roots only in a PREDECLARED local window."""
    scale=g.period/.3
    x=np.linspace(.65*scale,.74*scale,samples)
    values=np.array([jets(float(y),g)["p1"] for y in x])
    indexes=np.flatnonzero(values[:-1]*values[1:]<0)
    return [float(brentq(lambda z:jets(z,g)["p1"],float(x[j]),
                  float(x[j+1]),xtol=1e-14)) for j in indexes]


def clearance(g: Geometry, x: float, samples: int = 2001):
    """Sampled non-local clearance of the circular planar interface."""
    d=jets(x,g)
    R=d["w"]/(2*math.cos(d["apparent_angle_rad"]))
    yy=np.linspace(0,d["w"]/2,samples)
    if R*R<(d["w"]/2)**2:
        return -math.inf
    xx=x+R*math.sin(d["apparent_angle_rad"])-np.sqrt(
        np.maximum(0,R*R-yy*yy))
    k=2*math.pi/g.period
    uu=k*xx
    wall=g.nominal_width+g.a*np.cos(uu)+g.b*(
        np.cos(5*uu)-np.cos(7*uu))
    # The final point is at the prescribed wall contact: clearance zero.
    if abs(float(wall[-1]/2-yy[-1]))>1e-8:
        raise ValueError("The circular interface misses the wall endpoint")
    return float(np.min(wall[:-1]/2-yy[:-1]))


def exact_low_order_correlations(max_order: int = 5):
    """Exact integer-Fourier parity obstruction for every n-point offset."""
    if not (1<=max_order<=5):
        raise ValueError("Only orders <=5 are certified by the 5/7 pair")
    states={(0,0)}
    count=[]
    for n in range(1,max_order+1):
        states={(total+freq, (odd+int(abs(freq)>1))%2)
                for total,odd in states for freq in (-7,-5,-1,1,5,7)}
        if (0,1) in states:
            raise AssertionError(f"An odd-high-mode zero sum exists at n={n}")
        count.append(len(states))
    # It MUST fail at sixth order: (1,1,1,1,1,-5).
    states={(total+freq, (odd+int(abs(freq)>1))%2)
            for total,odd in states for freq in (-7,-5,-1,1,5,7)}
    assert (0,1) in states
    return count


def width_separation_um(b: float):
    u=np.linspace(0,2*math.pi,100001)
    return float(1000*2*abs(b)*np.max(abs(np.cos(5*u)-np.cos(7*u))))


def pair(g0: Geometry, b_minus: float, b_plus: float):
    """Return paired local onsets or explicit root-count failure."""
    gm=Geometry(g0.a,g0.period,g0.theta_deg,-b_minus,g0.nominal_width)
    gp=Geometry(g0.a,g0.period,g0.theta_deg,+b_plus,g0.nominal_width)
    rm,rp=stationary_roots(gm),stationary_roots(gp)
    if len(rm)!=1 or len(rp)!=1:
        return dict(ok=False,roots_minus=len(rm),roots_plus=len(rp))
    dm,dp=jets(rm[0],gm),jets(rp[0],gp)
    ok=(dm["p2"]<0 and dp["p2"]<0 and dm["B4"]<0<dp["B4"]
        and dm["jacobi_amplitude"]>0 and dp["jacobi_amplitude"]>0)
    return dict(ok=ok, roots_minus=1, roots_plus=1,
                s_minus_mm=rm[0],s_plus_mm=rp[0],
                B4_minus=dm["B4"],B4_plus=dp["B4"],
                p2_minus=dm["p2"],p2_plus=dp["p2"],
                jacobi_minus=dm["jacobi_amplitude"],
                jacobi_plus=dp["jacobi_amplitude"],
                clearance_minus_mm=clearance(gm,rm[0]),
                clearance_plus_mm=clearance(gp,rp[0]))


def stress_grid(b0: float = .00017):
    """405 discrete combinations, NOT a certified continuum enclosure."""
    rows=[]
    for a in (.049,.05,.051):
        for period in (.297,.3,.303):
            for angle in (30,31,32,33,34):
                g=Geometry(a,period,angle)
                for bm in (.8*b0,b0,1.2*b0):
                    for bp in (.8*b0,b0,1.2*b0):
                        outcome=pair(g,bm,bp)
                        rows.append(dict(a_mm=a,period_mm=period,
                            theta_deg=angle,b_minus_mm=bm,b_plus_mm=bp,
                            **outcome))
    return rows


def stationary_fold_creation(a=.049,period=.303,angle=34.):
    """Solve p'=p''=0 for positive-phase extra-root creation (local)."""
    def fun(q):
        g=Geometry(a,period,angle,float(q[1]))
        d=jets(float(q[0]),g)
        return [d["p1"],d["p2"]]
    answer=nonlinear_root(fun,[.70283,.000242],tol=1e-11)
    if (not answer.success or np.linalg.norm(answer.fun)>1e-7 or
            not (.0<answer.x[1]<.001)):
        raise ValueError("Fold-creation calculation failed")
    return dict(x_mm=float(answer.x[0]),b_mm=float(answer.x[1]),
                residual_norm=float(np.linalg.norm(answer.fun)))


def summary(b0=.00017):
    exact_low_order_correlations()
    nom=pair(Geometry(),b0,b0)
    if not nom["ok"]: raise ValueError("Nominal local onset class failed")
    cases=stress_grid(b0)
    if len(cases)!=405 or not all(r["ok"] and
        r["clearance_minus_mm"]>0 and r["clearance_plus_mm"]>0
        for r in cases):
        raise ValueError("Some grid points violate the local design gate")
    fold=stationary_fold_creation()
    return dict(model_scope="planar local + exact neutral-plate cylinders ONLY",
      correlation_orders_exact=list(range(1,6)),
      sixth_order_can_distinguish=True,
      nominal_b_um=b0*1000,nominal_theta_deg=32,
      max_pair_width_difference_um=width_separation_um(b0),
      nominal=nom,stress_points=len(cases),
      min_B4_margin=min(min(-r["B4_minus"],r["B4_plus"]) for r in cases),
      min_negative_p2_margin=min(min(-r["p2_minus"],-r["p2_plus"]) for r in cases),
      min_jacobi_amplitude=min(min(r["jacobi_minus"],r["jacobi_plus"]) for r in cases),
      min_sampled_clearance_mm=min(min(r["clearance_minus_mm"],
                                    r["clearance_plus_mm"]) for r in cases),
      worst_corner_fold_creation=fold,
      positive_amplitude_headroom_nm=(fold["b_mm"]-1.2*b0)*1e6)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-csv",action="store_true",
                       help="Print all 405 discrete stress-grid rows")
    args=parser.parse_args()
    if args.grid_csv:
        rows=stress_grid()
        writer=csv.DictWriter(sys.stdout,fieldnames=rows[0].keys())
        writer.writeheader();writer.writerows(rows)
    else:
        print(json.dumps(summary(),indent=2))


if __name__=="__main__":
    main()
