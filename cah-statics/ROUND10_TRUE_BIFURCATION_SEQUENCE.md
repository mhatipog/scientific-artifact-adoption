# Round 10 — A physically testable two-stage capillary bifurcation

**18 October 2026. Model-only research, not a fitted experiment.**

### Executive result

For a *fabrication-scale*, reflection-symmetric phase-paired channel

\[
w_\pm(x)=1+0.05\cos(kx)
   \pm 0.0004[\cos(2kx)-\cos(4kx)]\quad {\rm mm},
\qquad k=2\pi/(0.3\,{\rm mm}),
\]

at an assumed intrinsic sidewall Young angle 36.5°, planar
finite-volume Young–Laplace statics predicts:

| Observable | Negative Fourier phase A | Positive Fourier phase B |
|---|---:|---:|
| Symmetric local instability volume | 0.311037074 µL | 0.315944369 µL |
| Local onset | **Subcritical** | **Supercritical** |
| Nearby stable asymmetric branch after onset | None in examined branch | **Exists** |
| Subsequent asymmetric saddle-node | Not on this branch | **0.313349018 µL** |
| Centre displacement at second fold | — | **23.6007 µm** |
| Volume interval between continuous onset and fold | — | **2.59535 nL** |

The full 4-DOF static Hessian confirms the difference. This
suggests a *two-stage* response in B when volume decreases:
a continuous departure from symmetry, followed by loss of its
asymmetric local minimum. A instead loses its symmetric
local minimum at its subcritical onset.

**This is NOT a demonstrated movie of a jump.** The landing,
actual path, speed, and even realization of the predicted
instabilities depend on real capillary hysteresis, wetting
conditions, geometry, forcing and dynamics.

## 1. A previously unnoticed false-minimum mechanism

Earlier explorations eliminated cap shape by imposing the Young
angles at *both* contact lines at every displacement and then
eliminated plug half-length \(s=s(c)\) from fixed liquid volume:

\[
H(s-c)+H(s+c)=A_0.
\]

Along that **restricted** curve the exact energy derivative is

\[
\frac{dE_Y}{dc}
=\left[p(s-c)-p(s+c)\right]\,
\frac{2H'(s-c)H'(s+c)}
     {H'(s-c)+H'(s+c)}.
\]

**A zero of this derivative does NOT necessarily mean
\(p_L=p_R\).** It can also occur at a fold where
\(H'(s-c)=0\) or \(H'(s+c)=0\). Such a point is stationary
only in the restricted parametrization: varying the meniscus
cap angles independently at fixed total volume generally
unmasks the pressure imbalance.

For example, in phase A at its symmetric threshold,
the restricted curve has a false extremum near
\(c=33.416997\) µm with \(H'(s-c)\approx0\),
but
\[
\boxed{p_L-p_R\approx-0.07102\ {\rm mm}^{-1}\ne0.}
\]
It therefore FAILS the full 4-DOF capillary stationarity
conditions and must not be called a stable intermediate state.
The apparent restricted-curve extrema at even larger
displacements can also be of this kind. This error can
significantly mislead reported capillary-energy barriers.

**General methodological recommendation:** never infer stable
stick–slip branches from a one-coordinate constrained energy
curve alone when the elimination map has a singular
Jacobian. Test equal pressure, all cap-angle variations
and a projected full constrained Hessian.

## 2. Full equilibrium and second-variation test

Take the independent variables
\(z=(l,r,\theta_l,\theta_r)\), the exact planar area
\[
A=B(r)-B(l)+w(l)^2f(\theta_l)+w(r)^2f(\theta_r)
\]
and the normalized energy
\[
E=w(l)g(\theta_l)+w(r)g(\theta_r)
 -2\cos\theta_Y[S(r)-S(l)].
\]

Here \(B'=w\), \(S'=\sqrt{1+(w'/2)^2}\),
\(f(t)=(2t-\pi+\sin2t)/(8\cos^2t)\),
\(g(t)=(\pi-2t)/(2\cos t)\).

We independently solve
\[
p_L=p_R,\qquad A=A_0,
\]
using the natural Young angles only at actual stationary
states. We then verify
\(\|\nabla E+p\nabla A\|\) and diagonalize the full
4-variable capillary Hessian restricted to
\(\delta A=0\), which has **three** eigenvalues.
Both circular interfaces undergo a sampled wall-clearance
test; this is not a complete geometric/topological proof.

**Numerically checked states at \(\theta_Y=36.5^\circ\):**

- Phase A, just **1 nL above** symmetric onset:
  symmetric branch is strictly stable; a
  nonzero **index-one saddle** occurs at
  \(c=11.269843\) µm, with one negative and
  two positive constrained Hessian eigenvalues.
- Phase A, exactly at onset: symmetric branch
  has one zero eigenvalue; the first nearby asymmetric
  branches are unstable in the subcritical direction.
- Phase B, exactly at onset: a genuine **index-one
  saddle** occurs at \(c=33.629831\) µm.
  It is NOT the spurious \(H'=0\) extremum
  previously observed in the restricted calculation.
- Phase B, **1 nL below onset**: a *stable*
  asymmetric equilibrium exists at
  \(c\approx10.954\) µm, independently checked
  by three positive constrained eigenvalues.
- Phase B, secondary fold: the stable and
  index-one unstable asymmetric branches meet at
  \(c=23.600694\) µm, with the lowest Hessian
  eigenvalue zero to within the finite-difference
  tolerance.
- At the two onset volumes, both channels also have
  lower-energy **distant locally stable** equilibria
  at \(c=150\) µm. Their existence does NOT
  select which trajectory a dynamic experiment follows.

## 3. Energies relevant to measurable metastability

Use \(\gamma=0.072\) N/m and depth \(d=0.38\) mm
*only as illustrative conversion parameters*.
The physical surface energy scale is
\(\gamma d E_{\rm normalized}\).

- **Phase A**, at volume \(V_{A,*}+1\) nL,
  symmetric state to the verified index-one saddle:
  \(\Delta E=0.4289\) pJ.
- **Phase B**, at volume \(V_{B,*}\),
  symmetric critical state to the index-one saddle:
  \(\Delta E=6.9799\) pJ.

These differences are **candidate local saddle barriers**.
Identifying the *lowest global transition path* requires
following the saddle's downhill connections or a
minimum-energy-path algorithm. We have not established
global minimax barriers or thermally activated rates.
Contact-line friction, noise and hysteretic arrest
also remain unspecified.

The 7 pJ scale is far larger than molecular thermal
energy \(k_BT\) at room temperature. Therefore
pure thermal hopping between the deterministic wells
is not a plausible dominant experimental mechanism;
external perturbations, evaporative driving or
unmodeled defects are likely more important.

## 4. Model sensitivity to intrinsic wettability

The secondary fold was re-solved as a true
equal-pressure saddle-node for six assumed Young
angles, rather than by extrapolating one coefficient:

| Assumed θY | B continuous interval (nL) | B fold centre (µm) | (B fold − A onset) (nL) |
|---:|---:|---:|---:|
| 35.5° | 3.770 | 25.767 | 2.381 |
| 36.0° | 3.169 | 24.724 | 2.378 |
| 36.5° | 2.595 | 23.601 | 2.312 |
| 37.0° | 2.054 | 22.374 | 2.184 |
| 37.5° | 1.553 | 21.011 | 1.993 |
| 38.0° | 1.099 | 19.461 | 1.736 |

This difference is theoretically robust over this
angle sweep but the continuous interval narrows.
Resolving a ~1–4 nL evolution at total volumes
~310–318 nL needs careful volume calibration.

The **apparent in-plane angle** at the symmetric
critical right contact is only about
**11.1° in phase A** and **12.7° in phase B**
for the chosen intrinsic \(36.5^\circ\) angle.
This shallow contact geometry amplifies sensitivity
to fabrication imperfections and contact-line pinning.
It is **not** identical to the ~30° apparent
angle shown in the Poulose experiment.

## 5. Experimental falsification plan

The sharpest prediction is not whether a liquid
"jumps" at all. It is whether one phase-matched
channel exhibits a **resolvable continuous
symmetry-breaking segment spanning several nL**
before a secondary fold, while the other has
a subcritical onset under otherwise matched conditions.

For this controlled pair:

1. Independently measure fabricated width traces
   and full Fourier **phases**, not merely their
   common Fourier magnitudes; assess whether
   ±0.4 µm is resolved after manufacturing.
2. Independently calibrate **sidewall** intrinsic
   wettability, **plate** wettability and contact
   angle hysteresis under comparable conditions.
3. Record both contact-line positions and infer
   synchronized total liquid volume with
   uncertainties, then predict the two bifurcation
   sequences **blind**.
4. Include as baselines the original 2D sinusoidal
   theory and published competing fold/pitchfork
   predictions; measure whether phase information
   improves the out-of-sample result.
5. Validate with converged **nonneutral-plate 3D**
   contact-line eigenmodes and nonlinear quartic
   before asserting quantitative correspondence
   with a real rectangular PMMA chip.

## 6. Originality / prior art and readiness

Topography-induced folds, pitchforks and discontinuous
capillary motion are **not new phenomena**.
Relevant prior art includes
[Wells et al., Nature Communications 2018](https://doi.org/10.1038/s41467-018-03840-6)
and [Bisquert et al., Langmuir 2024](https://doi.org/10.1021/acs.langmuir.4c02908).
[Energy barriers in contact-angle hysteresis](https://www.sciencedirect.com/science/article/abs/pii/S000186860600131X)
have also been studied for decades.

What remains possibly distinctive is a **quantitative
comparison of phase-only wall perturbations at
identical roughness power spectrum**, showing
a different sequence of *verified free-cap-angle
full-Hessian bifurcations*. The new spurious-extremum
audit is an important safeguard against overclaiming
energy barriers.

**Professor verdict:**
- Internal planar theoretical verification: 7.5/10.
- Novelty confidence relative to the literature: 5/10.
- Experimental validation: 1/10.
- Submission readiness for *Langmuir* or
  *Physics of Fluids*: **4/10** pending measured
  tests and converged nonneutral 3D calculations.

## Reproducibility

Run:
```bash
pip install -r cah-statics/requirements.txt
python cah-statics/round10_full_equilibria.py
```
The script checks full constrained Hessians, normal
equilibria, real saddle-node conditions, a known false
restricted extremum, angle sensitivity and sampled
cap/wall clearance. It writes JSON/CSV and PNG/SVG
model-only outputs. Existing lightweight CAH PR CI
runs it; the heavy JAX 3D experiment remains manual.
