"""Exact high-order statistical degeneracy of local CAH onset (theory only).

Two reflection-symmetric channel widths are identical in their
translation-invariant correlations through FIFTH order, but have
opposite 2D nonlinear instability classes:
w_±(x)=1 + a*cos(k*x) ± b*[cos(5*k*x)-cos(7*k*x)].

This is not experimental validation or a nonlinear nonneutral-plate 3D result.
"""
import itertools, csv, json
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import numpy as np
import mpmath as mp
from scipy.optimize import brentq

A=.05; B=.000075; LAM=.3; THETA=36.
K=2*np.pi/LAM
OUTPUT=Path(__file__).parent/"outputs"
OUTPUT.mkdir(exist_ok=True)


def exact_invariants():
    freqs=[-7,-5,-1,1,5,7]
    results=[]
    for degree in range(2,7):
        tuples=[group for group in itertools.product(freqs,repeat=degree)
                if sum(group)==0]
        flip=[group for group in tuples
              if sum(abs(x) in (5,7) for x in group)%2]
        results.append(dict(degree=degree,
                            nonzero_frequency_product_count=len(tuples),
                            sign_flipping_count=len(flip),
                            sign_flipping_example=flip[0] if flip else None))
        assert bool(flip)==(degree==6)
    # Moment polynomials: constant coefficient of q(x)^degree,
    # q=a cos(kx)+b[cos(5kx)-cos(7kx)].
    modes={(1,1,0):Fraction(1,2),(-1,1,0):Fraction(1,2),
           (5,0,1):Fraction(1,2),(-5,0,1):Fraction(1,2),
           (7,0,1):Fraction(-1,2),(-7,0,1):Fraction(-1,2)}
    poly={(0,0,0):Fraction(1)}
    polys=[]
    for degree in range(1,7):
        new=defaultdict(Fraction)
        for (f,pa,pb),c in poly.items():
            for (g,qa,qb),q in modes.items():
                new[(f+g,pa+qa,pb+qb)]+=c*q
        poly=dict(new)
        p={(pa,pb):c for (f,pa,pb),c in poly.items() if f==0 and c}
        assert all(pb%2==0 for pa,pb in p) if degree<=5 else any(pb%2 for pa,pb in p)
        polys.append(p)
    assert polys[2]=={}
    assert polys[1]=={(2,0):Fraction(1,2),(0,2):Fraction(1)}
    difference=2*sum(float(c)*A**pa*B**pb for (pa,pb),c in polys[5].items() if pb%2)
    expected=3*A*B*(A**4-10*A*A*B*B+10*B**4)/8
    assert abs(difference-expected)<1e-22
    return results,difference


def critical_quartic(sgn,angle):
    """Independent pressure-root isolation and high-precision B4 formula."""
    b=sgn*B;theta=np.deg2rad(angle)
    def dp(x):
        w=1+A*np.cos(K*x)+b*(np.cos(5*K*x)-np.cos(7*K*x))
        wd=-K*(A*np.sin(K*x)+b*(5*np.sin(5*K*x)-7*np.sin(7*K*x)))
        wdd=-K*K*(A*np.cos(K*x)+b*(25*np.cos(5*K*x)-49*np.cos(7*K*x)))
        t=theta+np.arctan(wd/2)
        return 2*(-np.sin(t)*wdd/2/(1+(wd/2)**2)*w-np.cos(t)*wd)/w**2
    grid=np.linspace(.65,.74,2501)
    yy=dp(grid)
    crossings=np.where(yy[:-1]*yy[1:]<0)[0]
    assert len(crossings)==1, (sgn,angle,len(crossings))
    i=crossings[0]
    s=brentq(dp,grid[i],grid[i+1],xtol=2e-14)
    mp.mp.dps=32
    k=2*mp.pi/mp.mpf(str(LAM));a=mp.mpf(str(A))
    bb=mp.mpf(str(b)); th=mp.mpf(str(angle))*mp.pi/180
    sc=mp.mpf(str(s))
    def w(x):return 1+a*mp.cos(k*x)+bb*(mp.cos(5*k*x)-mp.cos(7*k*x))
    def wd(x):return -k*(a*mp.sin(k*x)+bb*(5*mp.sin(5*k*x)-7*mp.sin(7*k*x)))
    def t(x):return th+mp.atan(wd(x)/2)
    def f(t):return (2*t-mp.pi+mp.sin(2*t))/(8*mp.cos(t)**2)
    def H(x):
        return x+a*mp.sin(k*x)/k+bb*(mp.sin(5*k*x)/(5*k)-mp.sin(7*k*x)/(7*k))+w(x)**2*f(t(x))
    def p(x):return 2*mp.cos(t(x))/w(x)
    b4=(mp.diff(H,sc,2)*mp.diff(p,sc,2)-mp.diff(H,sc)*mp.diff(p,sc,3)/3)/4
    return float(s),float(b4)


def main():
    polyspectrum,difference=exact_invariants()
    # Bound: |cos(5u)-cos(7u)| = 2|sin(6u)sin(u)|
    # <=12 sin^2(u), so b<a/24 guarantees exact peak/trough 1±a.
    assert B<A/24
    x=np.arange(131072)*2*np.pi/131072
    v=[1+A*np.cos(x)+sgn*B*(np.cos(5*x)-np.cos(7*x))
       for sgn in (-1,1)]
    for stat in (np.mean,np.min,np.max,np.std):
        assert abs(stat(v[0])-stat(v[1]))<1e-12
    assert np.isclose(v[0].max(),1.05,atol=1e-12)
    assert np.isclose(v[0].min(),.95,atol=1e-12)
    np.testing.assert_allclose(np.abs(np.fft.rfft(v[0]-1))**2,
                               np.abs(np.fft.rfft(v[1]-1))**2,
                               rtol=1e-11,atol=1e-13)
    out=[]
    for sgn in (-1,1):
        s,b4=critical_quartic(sgn,THETA)
        out.append(dict(phase=sgn,critical_half_length_um=1000*s,
                        B4_local_2d=b4,
                        onset="subcritical" if b4<0 else "supercritical"))
    assert out[0]["B4_local_2d"]< -1400
    assert out[1]["B4_local_2d"]>700
    lo=brentq(lambda t:critical_quartic(-1,t)[1],30,31,xtol=1e-7)
    hi=brentq(lambda t:critical_quartic(1,t)[1],37,38,xtol=1e-7)
    assert 30.30<lo<30.32 and 37.68<hi<37.70
    maxdelta_um=1000*np.max(abs(v[1]-v[0]))
    details=dict(interpretation="2D analytical counterexample, no experimental or nonneutral 3D validation",
                 wall_width_mm="1 + .05*cos(kx) ± .000075*(cos(5kx)-cos(7kx))",
                 same=["mean","min/max","RMS","entire power spectrum",
                       "ALL translation-invariant correlations of order <=5"],
                 maximum_profile_difference_um=float(maxdelta_um),
                 first_differing_centered_moment_order=6,
                 sixth_centered_moment_difference_mm6=difference,
                 angle_window_deg=[lo,hi],quartic_at_36_degrees=out,
                 polyspectrum_invariance=polyspectrum)
    (OUTPUT/"highorder_isospectral.json").write_text(json.dumps(details,indent=2))
    with (OUTPUT/"highorder_isospectral.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(out[0]));wr.writeheader();wr.writerows(out)
    print("PASS exact two-through-five point spectral correlation equivalence")
    print("PASS first differing 6th moment",difference)
    print("PASS different B4 signs",out)
    print("PASS opposite-onset theta window",lo,hi)
    print("MODEL ONLY; no direct 3D B4/experiment comparison")

if __name__=="__main__":
    main()
