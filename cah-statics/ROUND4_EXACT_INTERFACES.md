# Round 4 — Exact stationary interfaces and physically smoothed walls

**Result status:** Independent Young–Laplace verification in two dimensions
and a provably exact *stationary* 3D benchmark. **No unrestricted 3D
stability calculation or experimental jump prediction has been completed.**

## Why a second solver matters

We previously called the planar circular-cap representation a "local
approximation". That wording was too broad. Under a constant-pressure,
constant-curvature, zero-gravity *two-dimensional* Young–Laplace equation,
the interface is precisely an arc of a circle. If the wall contact angles
are satisfied and the arc does not touch the solid elsewhere, it solves the
planar problem exactly. **The approximation is the reduction of an actual
finite-depth and potentially 3D meniscus to a planar cross-section**, not
the circularity of a genuine 2D constant-curvature solution.

For the two isospectral corrugations of Round 3, we separately solved

\[
x''(y)=\kappa[1+x'(y)^2]^{3/2}
\]

with the midplane symmetry condition \(x'(0)=0\), the actual sidewall
contact position \(x(w(s)/2)=s\), and sidewall slope
\(x'(w(s)/2)=\cot(\theta_Y+\alpha(s))\). Pressure \(\kappa\) is an
*unknown solved by the boundary-value problem*, not prescribed as an
answer. Initial curves and curvature were deliberately perturbed.

Both solutions agree with the circular-arc predictions to numerical
precision (max interface-position discrepancy < 1e-8 mm and curvature
discrepancy < 1e-9 mm^-1) and do not cross the channel walls.

**This corroborates the exact 2D stationary geometry. It does not
independently validate the bifurcation quartic coefficient, nor exclude
other 2D/3D modes.**

## What finite depth can be stated rigorously

If the upper and lower plates have an intrinsic contact angle of
**exactly 90°** and the sidewalls are extruded uniformly through depth
\(d\), a 2D circular cap lifted into an invariant cylinder satisfies the
full **3D stationary Young–Laplace equation**:

\[
\kappa_1=1/R,\qquad \kappa_2=0,\qquad
\Delta P=\gamma(\kappa_1+\kappa_2)=\gamma/R.
\]

The surface normal has zero vertical component and thus meets both
horizontal plates at 90°. Its sidewall contact angle is unchanged.
For the *z-invariant subspace*, volume and energy simply scale with
depth, so the relative sign of the quartic stability coefficient does
not depend on depth.

**Crucial limitation:** It has NOT been shown that these 3D cylinders
are stable to arbitrary depth-dependent perturbations. Non-neutral
(θ_top/bottom ≠ 90°) plates require a true 3D meniscus solution.
This is *not* a quantitative 3D validation of the Poulose channel.

## A quantitative fabrication/measurement diagnostic

For the isospectral test pair, physically Gaussian-smooth the **wall
profile** by a length \(\sigma\), attenuating each Fourier harmonic of
wavenumber \(nk\) by \(e^{-(nk\sigma)^2/2}\). This operation modifies the
actual walls, rather than just blurring their photographed image.

For assumed intrinsic \(\theta_Y=36.5^\circ\), the pair retains
opposite signs of \(B_4\) until

\[
\boxed{\sigma_* = 19.97694\;\mu\mathrm{m}.}
\]

After that, the positive-phase member changes sign and both are
subcritical (for the parameter range examined). This is a *model
manufacturing-sensitivity prediction*, not an experimental tolerance
measurement, and Gaussian smoothing is not a complete fabrication
process model. Image blur would not physically alter the bifurcation;
it would impair inference of the underlying phase.

## Scientific judgement after Round 4

- **Independent planar YL verification:** strong, but a known elementary
  consequence of uniform 2D curvature; not a standalone discovery.
- **Finite-depth result:** exact *stationary* solution for neutral
  horizontal plates, but stability and non-neutral plates unresolved.
- **Wall-geometry prediction:** quantitative response to deliberate
  physical smoothing. Potentially useful for microfabrication and
  metrology, but still untested against manufactured profiles.
- **Closest prior work:** Wells et al. (2018), Nature Communications,
  https://doi.org/10.1038/s41467-018-03840-6 and Bisquert et al. (2024),
  Langmuir, https://doi.org/10.1021/acs.langmuir.4c02908 already
  establish topography-driven snap bifurcations and their competition.

**Professor assessment:** 4/10 for journal submission (unchanged).
We have improved physical validity, not yet experimental strength or
demonstrated priority over closely related theory. A *Langmuir* or
*Physics of Fluids* paper could become credible after unrestricted
3D stability/finite-depth validation and an independently measured
sequence of jumps or deliberately manufactured paired channels.

**Next strict gate:** solve a full 3D capillary-surface stability problem
with freely moving sidewall contact lines, finite depth and realistic
top/bottom contact angles. Then compare predictions *blindly* with
independently measured wall topography, wetting and video trajectories.

Reproduce with:

```bash
python cah-statics/young_laplace_depth.py
```

GitHub Actions (PR/manual only) checks the BVP agreement and smoothing
threshold, and archives the model-only PNG/SVG and CSV/JSON data.
