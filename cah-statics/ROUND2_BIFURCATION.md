# Round 2 — Symmetry breaking distinguishes smooth shifts from jumps

**Status: new analytical/reproducible reduced-model result; NO quantitative experimental validation.**

## Reduction to one curvature and one volume primitive

For symmetric, sinusoidally corrugated walls, write the left and right front
locations as \(c-s\) and \(c+s\), where \(c\) is the plug centre, \(s\) its
half-length. With \(w(x)=w_0+a\cos(kx)\), constant intrinsic Young angle
\(\theta_Y\) and local wall inclination
\(\alpha(x)=\arctan[w'(x)/2]\), define

\[
p(x)=\frac{2\cos[\theta_Y+\alpha(x)]}{w(x)},\qquad
H(x)=\int_0^x w(t)\,dt+w(x)^2 f[\theta_Y+\alpha(x)]
\]

with \(f(\theta)=(2\theta-\pi+\sin2\theta)/(8\cos^2\theta)\).
For the two planar circular menisci (not a full 3-D Young–Laplace solution),

\[
V(c,s)/d=H(s-c)+H(s+c),\quad
p(s-c)=p(s+c).
\]

The full constrained-energy Hessian, including interface angles, has been
checked independently against the simple condition \(p'(s_*)=0\) on this
selected symmetric branch. At the critical volume, the fixed-volume energy
expansion is

\[
E(c;V_*)=E_0+B_4c^4+O(c^6), \qquad
\boxed{B_4=\tfrac14[H''p''-\tfrac13H'p''']_{s_*}}.
\]

The sign of \(B_4\) distinguishes **supercritical** (\(B_4>0\), stable
asymmetric branches appear continuously) from **subcritical** (\(B_4<0\),
nearby asymmetric branches are unstable, finite jumps become possible).
This is a *local* bifurcation criterion. The model does not give jump speeds
or unique landing selection without a dynamical/noise/defect assumption.

### Reproducible numerical results

Using \(w_0=1\) mm, ripple wavelength 0.30 mm and peak-to-peak width
modulation \(2a/w_0=10\%\), the quartic coefficient crosses zero at the
**assumed intrinsic angle 37.145146°**. For 30°, it is positive; for
40°, negative. The corresponding small asymmetric branches from the
full equal-pressure, equal-volume equations have positive versus negative
constrained stiffness respectively. Analytical \(B_4\) agrees with
an independent **45-digit, elliptic-integral energy** calculation.

| Peak-to-peak width modulation | Young angle at \(B_4=0\) |
|---|---:|
| 5% | 20.54585° |
| 7% | 27.29409° |
| 10% | 37.14515° |
| 14% | 49.92630° |
| 20% | 69.82103° |

The relevant critical positions lie at ~0.69–0.70 mm half-length for
the selected phase. The sign boundary is a property of local periodic
geometry and does not depend on the number of **complete periods** separating
the menisci, within this ideal infinite-channel model. That invariance
must be tested in finite experimental chips.

**Crucial correction to our earlier claim:** an unstable symmetric
equilibrium does *not always* imply a finite contact-line jump. At 30°
the local symmetry breaking is continuous. The earlier 209.7 µm jump
is a selected equal-volume landing at 40°; it is **not a universal
prediction for the published 30° apparent-angle images**.

### Experimental relevance and principal caveats

Poulose et al. (2023) report corrugation and intermittent motion, plus
an apparent angle near 30°. Their apparent angle is **not an independent
measurement of intrinsic Young angle**. Their wall trace is not available
here as numerical coordinates. Thus 37.145° is **not** a calibrated
transition angle for their experiment. Wall roughness, chemical pinning,
top/bottom 3D curvature and the validity of local circular arcs can shift
or remove the predicted symmetry-breaking transition.

The Yan et al. (2025) pore experiment provides higher-resolution pressure
and imagery, but uses a different throat–pore geometry: these numerical
thresholds must **not** be transplanted into its pressure curves.

**Next validation gate:** extract independent wall profile and intrinsic
wetting properties; solve the actual 2D Young–Laplace boundary value
problem (not local arcs), then blind-predict both interfaces event-by-event
with full uncertainty budgets. Plot measured \(\Delta c\) versus
\(V-V_c\) and check whether one finds \(\Delta c\propto|V-V_c|^{1/2}\)
(continuous onset) or a finite jump/hysteresis (discontinuous onset).


## Critical novelty audit — closest previous work

Two closely related studies **must** be acknowledged before treating these
bifurcations as novel:

1. **Wells, Ruiz-Gutiérrez, Le Lirzin, Nourry, Orme, Pradas &
   Ledesma-Aguilar, Nature Communications 9:1380 (2018)**,
   [Snap evaporation of droplets on smooth topographies](https://doi.org/10.1038/s41467-018-03840-6).
   This work already predicts and observes lateral snap transitions on smooth
   corrugated surfaces via a **subcritical pitchfork that precedes a
   symmetric fold**, including symmetric and asymmetric branches.
2. **Bisquert, Ruiz-Gutiérrez, Pradas & Ledesma-Aguilar,
   Langmuir 40:24387–24396 (2024)**,
   [Competing Bifurcations Determine Symmetry Breaking During Droplet
   Snaps on Smooth Patterned Surfaces](https://doi.org/10.1021/acs.langmuir.4c02908).
   This work explicitly compares fold and pitchfork competitions across
   topography amplitude, equilibrium contact angle, and growth/evaporation.

**Therefore: neither a pitchfork bifurcation, subcritical snap, nor the
concept of a topography/contact-angle parameter map is new.** Marc Pradas
coauthored both studies and the 2019 microchannel statics paper. Our
provisional incremental result is a particularly compact
\`H\`-and-\`p\` derivative criterion to distinguish *supercritical from
subcritical* onset for a different physical confinement (a finite plug
between **two symmetric corrugated walls** and two coupled free interfaces).

The numerical value 37.145° has no experimental standing without
independently measured wall profiles and Young angles. We do **not** claim
the coefficient is absent from the wider literature; a detailed
method-by-method comparison with the above authors is still required.

**Important mathematical qualification:** in *planar* capillarity an
unobstructed, constant-curvature circular arc satisfying the wall
contact-angle condition is an exact Young–Laplace solution within
its symmetric interface class. However, the current model has not
excluded additional asymmetric contact-line modes on opposite walls,
nonlocal contact of an arc with a corrugation, or the experimental
channel's fully three-dimensional meniscus response. Consequently the
full experimental stability threshold still needs validation.

**Professor-style publishability conclusion (17 October 2026):**
not journal-ready. As an internal reproducible mathematical note it is
useful; for *Langmuir* or *Physics of Fluids* it needs demonstrably new
predictive value over Wells (2018) and Bisquert (2024), preferably
blind quantitative comparison to measured corrugated-channel trajectories.
*Physical Review Fluids* would additionally need a broad or universal
insight beyond the special sinusoidal example. A repo update and
potential collaborator discussion are appropriate now; a submission is not.


### Scholarly positioning (critical)

Symmetry breaking, capillary metastability, and Haines jumps are known
topics, and there are contemporary quantitative experiments (Yan et al.,
JCIS 2025; Zhang et al., Water Resources Research 2026). The *potential*
original contribution is the compact local **quartic-sign criterion**
and its geometry–wettability map for two interacting free menisci,
not discovering jumps or multistability.

A stand-alone new journal article is **premature** until an exact-interface
comparison and at least one experimental sequence survive prediction with
independently measured parameters. A technical preprint or short theory
note might be justified after full literature novelty screening and
independent numerical boundary-value validation.

## Reproducibility

```sh
pip install -r cah-statics/requirements.txt
python cah-statics/bifurcation.py
```

The script regenerates `quartic_sign_change.csv`,
`pitchfork_exact_branches.csv`, and SVG/PNG phase-diagram figures
under `cah-statics/outputs/`. GitHub Actions executes it on relevant
pull requests and by manual dispatch, without a schedule.

Sources:
- [Hatipogullari et al., PRF 2019](https://doi.org/10.1103/PhysRevFluids.4.044008)
- [Poulose et al., 2023](https://arxiv.org/abs/2210.16962)
- [Yan et al., JCIS 2025](https://doi.org/10.1016/j.jcis.2025.137884)
- [Sun & Santamarina, PRE 2019](https://doi.org/10.1103/PhysRevE.100.023115)
- [Zhang et al., WRR 2026](https://doi.org/10.1029/2024WR039565)
