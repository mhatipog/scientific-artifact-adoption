"""Exact single-front-to-two-front local quartic bridge and measurement sign bound.

Model premise: Q=H(s) is strictly increasing, p=p(Q) has three derivatives,
p_Q(Qc)=0, and the symmetric plug has the fixed-volume coefficient
B4=(H'' p''-H' p'''/3)/4 as in the stated reduced model.
This is algebra, not a 3D stability or experiment claim.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Bridge:
    pressure_curvature_in_Q: float
    pressure_third_derivative_in_Q: float
    quartic_from_volume: float
    quartic_from_geometry: float
    identity_residual: float


def bridge_at_stationary(H1: float, H2: float, p1: float,
                         p2: float, p3: float,
                         pressure_stationarity_tolerance: float = 1e-8) -> Bridge:
    """Apply the exact chain rule at p'(s)=0, H1=H'(s)>0.

    These formulas for p_QQ and p_QQQ use stationarity and are *not*
    valid off the critical point without additional p1 and H3 terms.
    """
    values=(H1,H2,p1,p2,p3)
    if not all(math.isfinite(x) for x in values) or H1<=0:
        raise ValueError("H1 must be positive; all jets finite")
    if (pressure_stationarity_tolerance<0 or
        abs(p1)>pressure_stationarity_tolerance):
        raise ValueError("Pressure is not stationary within declared tolerance")
    pQQ=p2/H1**2
    pQQQ=(p3-3*(H2/H1)*p2)/H1**3
    b4_vol=-H1**4*pQQQ/12
    b4_geom=(H2*p2-H1*p3/3)/4
    return Bridge(pQQ,pQQQ,b4_vol,b4_geom,b4_vol-b4_geom)


@dataclass(frozen=True)
class SignCertificate:
    third_derivative_estimate: float
    worst_case_error_bound: float
    third_derivative_sign: str
    local_quartic_sign: str


def pressure_third_derivative_sign_bound(
    p_minus_2h: float, p_minus_h: float,
    p_plus_h: float, p_plus_2h: float,
    h: float, pressure_absolute_error: float,
    certified_fifth_derivative_bound: float
) -> SignCertificate:
    """Conditional, rigorous four-pressure-point derivative sign test.

    Suppose each measured pressure has absolute error <= epsilon; p(Q)
    is C^5 on [Qc-2h,Qc+2h], and |d^5p/dQ^5| <= M5 there. Then

      |D3 - p'''(Qc)| <= 3*epsilon/h^3 + (17/60)*M5*h^2.

    The last constant follows by taking absolute fifth-order Taylor
    remainders individually at ±h and ±2h; no probabilistic independence
    or stochastic noise assumption is used. This is a sufficient bound,
    deliberately not optimized. It does not certify p_Q(Qc)=0, which
    must be separately established, nor the reduced-model assumptions.
    """
    data=(p_minus_2h,p_minus_h,p_plus_h,p_plus_2h,
          h,pressure_absolute_error,certified_fifth_derivative_bound)
    if (not all(math.isfinite(x) for x in data) or h<=0 or
        pressure_absolute_error<0 or certified_fifth_derivative_bound<0):
        raise ValueError("Nonfinite inputs or invalid positive/error bounds")
    derivative=(p_plus_2h-2*p_plus_h+2*p_minus_h-p_minus_2h)/(2*h**3)
    bound=(3*pressure_absolute_error/h**3+
           (17/60)*certified_fifth_derivative_bound*h*h)
    if derivative>bound:
        p_sign="positive"
        b_sign="negative"
    elif derivative < -bound:
        p_sign="negative"
        b_sign="positive"
    else:
        p_sign="indeterminate"
        b_sign="indeterminate"
    return SignCertificate(derivative,bound,p_sign,b_sign)
