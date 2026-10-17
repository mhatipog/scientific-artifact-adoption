"""Analytical two-front symmetry-breaking test (no fit to experimental data).

Theory for a symmetric sinusoidal channel with local circular menisci:

 V/d = H(s+c) + H(s-c),  p(x) = 2 cos(theta + atan(w'/2))/w
 H(x) = integral_0^x w(t)dt + w(x)^2 f(theta + atan(w'/2)).

At p'(s*)=0 the symmetric translation mode softens.  At fixed V:
 E(c;V*) = E(0;V*) + B4*c^4 + O(c^6)
 B4 = [H''(s*) p''(s*) - H'(s*) p'''(s*) / 3] / 4.

B4>0: nearby asymmetric states stable on the side V<V* (supercritical).
B4<0: nearby asymmetric states unstable on the side V>V* (subcritical).
These statements are local and conditional on all remaining modes stable.

All results are MODEL predictions; no real channel image / jump is fitted.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.optimize import brentq

from model import PlugModel

OUT = Path(__file__).parent / "outputs"
OUT.mkdir(exist_ok=True)
mp.mp.dps = 40
PERIOD = mp.mpf(".3")
WIDTH = mp.mpf("1")


def coefficient(angle_deg, amplitude_mm):
    """Return critical half-length s*, quartic coefficient B4 and functions."""
    angle = mp.radians(angle_deg)
    a = mp.mpf(str(amplitude_mm))
    k = 2 * mp.pi / PERIOD

    def w(x):
        return WIDTH + a * mp.cos(k * x)

    def alpha(x):
        return mp.atan(-a * k * mp.sin(k * x) / 2)

    def cap_f(t):
        return (2 * t - mp.pi + mp.sin(2 * t)) / (8 * mp.cos(t)**2)

    def primitive(x):
        return WIDTH * x + a * mp.sin(k * x)/k

    def h(x):
        return primitive(x) + w(x)**2 * cap_f(angle + alpha(x))

    def pressure(x):
        return 2 * mp.cos(angle + alpha(x)) / w(x)

    # This tracks the selected critical point (2+fraction) periods from zero.
    s = mp.findroot(lambda x: mp.diff(pressure, x),
                    (mp.mpf("2.25") * PERIOD, mp.mpf("2.43") * PERIOD))
    if not mp.mpf("2.2") * PERIOD < s < mp.mpf("2.45") * PERIOD:
        raise AssertionError("Unexpected curvature extremum; check branch")
    hp = mp.diff(h, s)
    hpp = mp.diff(h, s, 2)
    p2 = mp.diff(pressure, s, 2)
    p3 = mp.diff(pressure, s, 3)
    quartic = (hpp * p2 - hp * p3 / 3) / 4
    return float(s), float(quartic), (h, pressure)


def fixed_volume_quartic_check(angle_deg=40., amplitude_mm=.05):
    """Independent high-precision check via exact elliptic wall arclength."""
    s0, b4, (h, p) = coefficient(angle_deg, amplitude_mm)
    s0 = mp.mpf(str(s0))
    a = mp.mpf(str(amplitude_mm))
    theta = mp.radians(angle_deg)
    k = 2 * mp.pi / PERIOD
    b = a*k/2
    def w(x): return WIDTH + a*mp.cos(k*x)
    def t(x): return theta+mp.atan(-a*k*mp.sin(k*x)/2)
    def arc(t_): return (mp.pi-2*t_)/(2*mp.cos(t_))
    def S(x): return mp.ellipe(k*x,-b*b)/k
    v0 = 2*h(s0)
    def E(c,s):
        x1,x2=s-c,s+c
        return w(x1)*arc(t(x1))+w(x2)*arc(t(x2))-2*mp.cos(theta)*(S(x2)-S(-x1))
    e0=E(mp.mpf("0"),s0)
    cvals=[mp.mpf(z) for z in [".008",".012",".016",".020",".024",".028"]]
    M=mp.matrix([[c**4,c**6,c**8,c**10,c**12,c**14] for c in cvals])
    differences=[]
    for c in cvals:
        root=mp.findroot(lambda s:h(s-c)+h(s+c)-v0,
                         (s0-mp.mpf(".002"),s0+mp.mpf(".002")))
        differences.append(E(c,root)-e0)
    fit=mp.lu_solve(M,mp.matrix(differences))
    return float(fit[0]),float(b4)


def near_branches(theta):
    """Numerically solve *both* interface pressures and the volume constraint."""
    m=PlugModel(theta_deg=theta)
    sc, _, _=coefficient(theta,.05)
    vc=float(m.volume(-sc,sc))
    rows=[]
    for dv in np.r_[np.linspace(-.0004,-.00001,25),
                    np.linspace(.00001,.0004,25)]:
        volume=vc+float(dv)
        def fixed_s(c):
            return brentq(lambda s:float(m.volume(c-s,c+s)-volume),
                          sc-.03,sc+.03,xtol=1e-13)
        s0=fixed_s(0.)
        state=np.array([-s0,s0,
              float(m.theta(-s0,"left")),float(m.theta(s0,"right"))])
        stable=min(m.constrained_stability(state)[0])>0
        rows.append((theta,float(dv)*1000,0.,"symmetric",stable))
        def residual(c):
            s=fixed_s(c)
            return float((m.p(c-s,"left")-m.p(c+s,"right"))/c)
        cs=np.linspace(.0004,.018,42)
        values=[residual(c) for c in cs]
        signs=np.where(np.array(values[:-1])*np.array(values[1:])<0)[0]
        if len(signs):
            ix=signs[0]
            c=brentq(residual,cs[ix],cs[ix+1],xtol=1e-13)
            s=fixed_s(c)
            z=np.array([c-s,c+s,float(m.theta(c-s,"left")),
                        float(m.theta(c+s,"right"))])
            ev,force_residual=m.constrained_stability(z)
            assert np.linalg.norm(force_residual)<1e-7
            assert abs(float(m.volume(c-s,c+s)-volume))<1e-10
            for sign in [-1,1]:
                rows.append((theta,float(dv)*1000,sign*c*1000,
                             "asymmetric",min(ev)>0))
    return rows


def run():
    # Normal-form coefficient independently validated from the energy itself.
    fitted, analytical = fixed_volume_quartic_check()
    assert abs((fitted-analytical)/analytical)<2e-6, (fitted, analytical)
    s30,b30,_=coefficient(30,.05)
    s40,b40,_=coefficient(40,.05)
    assert b30>0 and b40<0
    assert abs(b40 + 322.29545)<.015
    amp_brackets=[(.025,20,23),(.035,25,30),(.05,36,38),
                  (.07,48,51),(.1,68,72)]
    thresholds=[]
    for a,lo,hi in amp_brackets:
        cross=brentq(lambda t:coefficient(t,a)[1],lo,hi,xtol=1e-8)
        s,at_zero,_=coefficient(cross,a)
        assert abs(at_zero)<1e-4
        thresholds.append((200*a,cross,s))
    assert abs(thresholds[2][1]-37.1451464)<1e-5

    with (OUT/"quartic_sign_change.csv").open("w",newline="") as f:
        wr=csv.writer(f);wr.writerow(["width_corrugation_peak_to_peak_percent",
                                     "young_angle_at_B4_zero_deg","scritical_mm"])
        wr.writerows(thresholds)
    branches=near_branches(30)+near_branches(40)
    with (OUT/"pitchfork_exact_branches.csv").open("w",newline="") as f:
        wr=csv.writer(f);wr.writerow(["assumed_young_angle_deg","V_minus_Vcritical_nL",
                                     "plug_center_um","branch","stable_full_hessian"])
        wr.writerows(branches)
    # Check representative near-centre stability sign on both sides.
    assert any(row[0]==30 and row[3]=="asymmetric" and row[4] for row in branches)
    assert any(row[0]==40 and row[3]=="asymmetric" and not row[4] for row in branches)

    fig=plt.figure(figsize=(14.7,9.5))
    grid=fig.add_gridspec(2,2,left=.09,right=.96,top=.88,bottom=.17,hspace=.59,wspace=.27)
    a=fig.add_subplot(grid[0,:])
    xs=np.array([r[0] for r in thresholds]);ys=np.array([r[1] for r in thresholds])
    a.fill_between(xs,ys,85,alpha=.16,label="B4 < 0: subcritical")
    a.fill_between(xs,0,ys,alpha=.12,label="B4 > 0: supercritical")
    a.plot(xs,ys,"o-",lw=2.3,label="B4 = 0",color="#303956")
    a.scatter([10],[thresholds[2][1]],s=135,facecolors="white",
              edgecolors="#303956",linewidth=2.2,zorder=5)
    a.set(xlabel="Width corrugation (peak-to-peak % of nominal 1 mm width)",
          ylabel="Assumed intrinsic Young angle (degrees)",
          xlim=(4,21),ylim=(12,85),
          title="A | Analytic geometry–wettability boundary")
    a.text(15.1,70,"B4 < 0  |  finite jumps possible",
           fontsize=13,weight="bold")
    a.text(10,17,"B4 > 0  |  smooth symmetry breaking",
           fontsize=13,weight="bold",ha="center")
    a.grid(alpha=.14);a.legend(frameon=False,loc="upper left")

    for j,ang in enumerate([30,40]):
        ax=fig.add_subplot(grid[1,j])
        for branch,stable,color in [("symmetric",True,"#15867d"),
                                    ("symmetric",False,"#bb7747"),
                                    ("asymmetric",True,"#15867d"),
                                    ("asymmetric",False,"#bb7747")]:
            rows=[r for r in branches if r[0]==ang and r[3]==branch and r[4]==stable]
            if rows:
                ax.scatter([r[1] for r in rows],[r[2] for r in rows],
                           s=17,color=color,
                           label=branch+" "+("stable" if stable else "unstable"))
        ax.axvline(0,linestyle=":",linewidth=1)
        ax.set(xlabel="V − Vcritical (nL)",
               ylabel="Plug centre (µm)" if j==0 else "",
               title=("B" if j==0 else "C")+f" | {ang}°: "+
                 ("nearby stable branches" if ang==30 else "nearby unstable branches"))
        ax.grid(alpha=.13);ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Two menisci, two distinct instability mechanisms",
                 y=.965,fontsize=22,weight="bold")
    fig.text(.5,.08,"THEORETICAL TEST · sinusoidal walls, local circular arcs, no measured trajectory or wettability fit",
             ha="center",fontsize=10)
    fig.savefig(OUT/"symmetry_breaking_phase_diagram.svg")
    fig.savefig(OUT/"symmetry_breaking_phase_diagram.png",dpi=175)
    plt.close(fig)
    print("PASS: analytic quartic coefficient matches independent fixed-volume energy to 2e-6")
    print("PASS: exact equilibria confirm stable/unstable near-onset branches")
    print("10% corrugation: B4=0 at Young angle",thresholds[2][1],
          "; B4(30°)=",b30,"; B4(40°)=",b40)
    print("MODEL PREDICTIONS ONLY: experimental agreement remains untested.")


if __name__ == "__main__":
    run()
