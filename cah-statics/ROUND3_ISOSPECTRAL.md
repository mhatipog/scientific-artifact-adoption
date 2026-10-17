# Round 3 — Isospectral wall shapes with opposite instability classes

**Computational finding, 17 October 2026. Theory only; no experimental matching.**

## A controlled geometry test

Replace the assumed sinusoidal channel by a pair of profiles

\[
w_{\pm}(x)=w_0+a\cos(kx)\pm b[\cos(2kx)-\cos(4kx)],
\]

with \(w_0=1\) mm, \(a=0.05\) mm, \(b=0.0005\) mm and
\(k=2\pi/(0.3\,\mathrm{mm})\). These are both *reflection-symmetric*
profiles, so the symmetric pitchfork theory applies.

The two channels have the **same** period, mean width, minimum and
maximum widths (0.95 and 1.05 mm), RMS width variation, RMS slope and
entire Fourier magnitude/power spectrum. The relative **Fourier phases**
are different. In particular, the two small 2nd/4th harmonic corrections
change sign simultaneously. This is not a coordinate translation, since
the fundamental remains fixed.

## Reproducible calculation at an ASSUMED θY = 36.5°

| Input/prediction | Minus harmonic phase | Plus harmonic phase |
|---|---:|---:|
| Geometry correction b | −0.5 µm | +0.5 µm |
| Critical meniscus half-length | 700.331 µm | 697.477 µm |
| Critical liquid volume | 0.310270 µL | 0.315995 µL |
| Fourth-order coefficient B4 | **−1080.33** | **+531.77** |
| Adjacent asymmetric branches | **unstable** | **stable** |
| Local transition | **subcritical** | **supercritical** |

For this isospectral pair, the signs stay opposite for intrinsic Young
angles approximately **33.22008° to 39.30682°**, a 6.09°-wide interval.
The two endpoint values are located by bracketing the zero of the same
quartic coefficient, separately for each shape; the included script
recomputes and asserts this window. This is still a **theoretical angle
window**, not an experimental error bar.

The coefficient is computed via the analytic criterion from Round 2,

\[
B_4=\frac14\left(H''p''-\frac13H'p'''\right)_{s_c},
\qquad p'(s_c)=0,
\]

not a fit. We additionally solve the *full equal-pressure/equal-volume*
conditions for nearby asymmetric states, evaluating constrained energy
Hessian eigenvalues. The branch stability agrees with the criterion.
Both circular interfaces clear the wall in the tested states.

**The result is much stronger than matching only roughness amplitudes**:
the power spectrum is rigorously identical, while the relative spatial
phase of the pattern determines the local bifurcation class.

### Why this is not a full claim of smooth versus jumping motion

Positive B4 establishes locally stable asymmetric states at bifurcation,
not the global absence of other discontinuities or distant metastable
branches. Negative B4 establishes locally unstable emerging branches
and *permits* a finite jump, but static theory alone does not select the
landing state or jump speed. We should say "opposite local onset classes"
rather than "one channel always jumps and the other never jumps."

### Why the 0.5 µm scale matters

The derivative criterion contains \(p'''\); higher spatial harmonics can
have a disproportionately large effect because differentiation weights
each wavenumber by powers of k. Roughness descriptors based on
\(\mathrm{RMS}(w)\), \(\mathrm{RMS}(w')\), or the power spectrum alone
can miss the **relative-phase** information needed to predict this
localized instability. The model therefore offers a *metrology*
hypothesis: channel-profile phase data may be indispensable for
deterministic meniscus-jump prediction.

### Caveats and strongest novelty competitors

Prior work **already** shows topography-dependent capillary bifurcations:
[Wells et al., Nature Communications (2018)](https://doi.org/10.1038/s41467-018-03840-6);
[Bisquert et al., Langmuir (2024)](https://doi.org/10.1021/acs.langmuir.4c02908).
Also relevant are
[Iliev et al., PRE (2018)](https://doi.org/10.1103/PhysRevE.97.042801)
and [Iliev et al., PRE (2023)](https://doi.org/10.1103/PhysRevE.107.024802),
which analyze wetting hysteresis and roughness factors. We have **not**
established that phase sensitivity in a full meniscus or in broader
surface statistics is new. Our particular isospectral counterexample
is a useful specific falsifiable prediction but not yet a journal
article.

This is a **plan-view 2D constant-curvature model**. Finite-depth effects
may change the branch structure, and nonlocal interactions or extra
contact-line modes may destabilize assumed equilibria. More
importantly, the channels used in Poulose et al. (2023) have NOT been
measured here as numerical wall profiles; the 36.5° is an arbitrary
intrinsic angle, not their apparent 30°.

## Recommended next result for publication

Use *real measured wall profiles* (ideally submicron accurate) to predict
the entire sequence of jump positions, with separately calibrated
wettability and no fitted jump lengths. Solve true 2D Young–Laplace
boundary-value problems, then extend to finite depth and compare
quantitatively. The best case is fabricating an **isospectral pair**
and measuring different snap onset classes directly.

**Professor assessment, round 3:** model insight 7/10, novelty confidence
4/10, experimental strength 1/10, submission readiness **4/10**.
Possible *Langmuir* or *Physics of Fluids* only after rigorous
independent validation; *PRF* requires stronger generality and impact.

## Reproduction

Run `python cah-statics/phase_isospectral.py` after installing the
requirements. GitHub Actions (manual or pull-request on changed paths)
regenerates `isospectral_phase_switch.csv` and SVG/PNG figures.
There is **no scheduled action or experimental data redistribution**.
