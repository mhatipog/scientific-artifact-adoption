"""A local non-identifiability theorem for capillary onset from linear response.

At a pressure-stationary circular meniscus, add an even C-infinity bump whose
jets w,w',w'',w''' vanish at the stationary position while w'''' changes.
For a symmetric two-meniscus channel, apply reflected bumps at both +/-s.
The exact neutral-cylinder Jacobi spectrum, neutral plate susceptibility,
critical position and lower pressure jets do not change, but the reduced
fixed-volume quartic B4 changes *linearly*. This is a model counterexample,
not a proposed nanometre-scale fabrication procedure.
"""
from __future__ import annotations
import argparse
import json
import math
from spectral_isolation import Channel, critical_state


def fourth_jet_sensitivity(state: dict, theta_deg: float) -> dict:
    """Exact dB4/d(w'''') in the planar circular-interface reduced model.

    The physical wall-curvature Robin data use only w,w',w''.
    p' uses w'', p'' uses w''', and p''' first uses w''''.
    """
    w,w1,w2=(state[x] for x in ("w","w1","w2"))
    theta=math.radians(theta_deg)
    z=w1/2
    den=1+z*z
    alpha=theta+math.atan(z)
    z1=w2/2
    a1=z1/den
    c=(2*alpha-math.pi)/8
    sec2=1/math.cos(alpha)**2
    hf=c*sec2+math.tan(alpha)/4
    hf1=sec2/2+2*c*sec2*math.tan(alpha)
    Hp=w+2*w*w1*hf+w*w*hf1*a1
    if Hp <= 0 or math.sin(alpha) <= 0:
        raise ValueError("This sign proof requires H'>0 and sin(apparent angle)>0")
    dp3_d_w4=-math.sin(alpha)/(w*den)
    dB4_d_w4=-Hp*dp3_d_w4/12
    return dict(H_prime_mm=Hp,theta_apparent_deg=math.degrees(alpha),
                dp3_d_w4=dp3_d_w4,dB4_d_w4=dB4_d_w4)


def invisible_flip(state: dict, theta_deg: float, factor: float = 2.0) -> dict:
    """Set delta w'''' so local B4 is (1-factor)*B4_initial."""
    if factor <= 1:
        raise ValueError("factor>1 required for a strict sign reversal")
    slope=fourth_jet_sensitivity(state,theta_deg)["dB4_d_w4"]
    baseline=state["B4"]
    delta=-factor*baseline/slope
    return dict(original_B4=baseline,delta_w4_per_mm3=delta,
                modified_B4=baseline+slope*delta,
                critical_position_mm=state["s_mm"],sensitivity=slope,
                neutral_cylinder_unchanged=True,
                full_3d_linear_spectrum_unchanged=True,
                plate_susceptibility_unchanged=True,
                caveat="Only exact for infinitesimally localized smooth jet modifications; global folds not asserted")


def symmetric_bump(x: float, s: float, radius: float,
                   delta_w4: float, period_mm: float | None = None) -> float:
    """Reflection-even smooth physical-width perturbation, not its derivative.

    The centered profile (delta_w4/24)*u^4*exp(1-1/(1-(u/r)^2))
    vanishes with every derivative at |u|=r, and has fourth derivative
    delta_w4 at u=0. Optional periodization uses only nearby centers.
    """
    if radius <= 0:
        raise ValueError("Bump radius must be positive")
    if period_mm is not None and period_mm <= 2*radius:
        raise ValueError("Periodic copies may overlap")
    values=0.0
    shifts=(0,) if period_mm is None else range(
        math.floor((x-abs(s)-radius)/period_mm)-1,
        math.ceil((x+abs(s)+radius)/period_mm)+2)
    for n in shifts:
        for center in (-s,s):
            u=x-center-n*(period_mm or 0)
            t=u/radius
            if abs(t)<1:
                values+=(delta_w4/24)*u**4*math.exp(1-1/(1-t*t))
    return values


def report() -> dict:
    g=Channel()
    states=[critical_state(phase,g) for phase in (-1,1)]
    checks=[]
    for state in states:
        example=invisible_flip(state,g.theta_deg,2)
        radius=.004
        s=state["s_mm"]
        samples=[s+radius*j/200 for j in range(-200,201)]
        excursion=max(abs(symmetric_bump(x,s,radius,
                     example["delta_w4_per_mm3"])) for x in samples)
        example["illustrative_bump_radius_mm"]=radius
        example["illustrative_max_width_change_nm"]=excursion*1e6
        checks.append(example)
    return dict(statement="same neutral-plate 3D linear operator and Q, opposite reduced quartic",
                geometry_pair=[{"phase":sgn,**r} for sgn,r in zip((-1,1),checks)],
                status="exact local jet theorem in the reduced two-front model; not global nonneutral 3D evidence")


if __name__=="__main__":
    print(json.dumps(report(),indent=2))
