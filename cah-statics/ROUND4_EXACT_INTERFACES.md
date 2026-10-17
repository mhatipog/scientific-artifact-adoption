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


## Exact 3D spectral reduction at neutral plates — next numerical gate

There is a sharper statement than equilibrium extrusion alone. For the
z-invariant cylindrical meniscus at 90° plates, its Jacobi (second
variation) operator is separable in planar interface arclength \(s\)
and the vertical coordinate \(z\):

\[
\mathcal J=-\partial_s^2-\partial_z^2-\kappa^2.
\]

The *moving contact lines* impose capillary Robin conditions at the
extruded sidewalls; their coefficients depend on actual wall curvature
and Young angle, but not on \(z\). Neutral horizontal plates impose
Neumann vertical conditions. Thus writing
\(u(s,z)=\psi_j(s)\cos(n\pi z/d)\) yields exactly

\[
\boxed{\lambda_{j,n}=\mu_j+(n\pi/d)^2},
\]

where \(\mu_j\) are eigenvalues of the **unconstrained** planar capillary
Jacobi operator *including its free-contact-line Robin conditions*.

A fundamental constraint distinction: the fixed *total* liquid-volume
constraint acts on the \(n=0\) family, but first-order volume variation
of every \(n\ge1\) vertical mode vanishes automatically. Consequently
**a stable volume-constrained planar branch is not sufficient by
itself to prove 3D stability**. The nonuniform modes are all stable
provided

\[
\mu_{\min}+(\pi/d)^2>0.
\]

If \(\mu_{\min}<0\), a sufficient depth threshold is
\(d<\pi/\sqrt{-\mu_{\min}}\); if \(\mu_{\min}\ge0\) the vertical modes
are positive at every depth. This applies only under the stated
idealized geometry/neutral-plate assumptions and is not an assertion
that the experimental chip satisfies those boundary conditions.

### A zero-mode theorem at the critical state

We can sharpen the previous *conditional* criterion **at the two
particular isospectral critical points**. Write the right cap as
\[
x(y;s)=s+R(s)\sin\vartheta(s)-\sqrt{R(s)^2-y^2},\quad
R(s)=\frac{w(s)}{2\cos\vartheta(s)},\quad
\vartheta(s)=\theta_Y+\arctan[w'(s)/2].
\]
The local pitchfork threshold is \(p'(s_*)=0\); since \(p=1/R\),
this means \(R'(s_*)=0\). Differentiating the *family of exact
Young–Laplace equilibria* at fixed interior \(y\) gives
\[
\left.\partial_s x(y;s)\right|_{s_*}
=A,\qquad
A=1+\frac{w(s_*)}{2}\,\vartheta'(s_*).
\]
Thus the normal deformation is
\[
u(y)=A\sqrt{1-(y/R)^2}.
\]
It solves the *homogeneous* capillary Jacobi equation, including the
moving-contact-line boundary conditions, since it comes from the
derivative of exact Young-angle equilibria at the pressure extremum.
If \(A\ne0\), this function has **no interior zeros**, so by the
Sturm oscillation theorem it is the *principal* planar Jacobi
eigenfunction and \(\mu_{\min}=0\) at that critical state.

For the two isospectral profiles at assumed \(36.5^\circ\),
\[
A_- = 3.4797654,\qquad A_+ = 2.9235464.
\]
Both are strictly positive, and an independent central
finite-difference of the entire family of circular interfaces
reproduces the constant derivative to better than \(2\times10^{-6}\)
mm per mm of contact-line displacement in the tested interior.

It follows that **at the selected critical state**, every
\(n\ge1\) depth-varying mode of the freely moving circular cap
satisfies
\[
\lambda_{j,n}\ge (\pi/d)^2>0.
\]
For \(d=0.38\) mm this lower bound is approximately
\(68.34\,\mathrm{mm}^{-2}\), with the Jacobi operator expressed in
physical arclength units. The identical conclusion applies to
the reflected left meniscus; the two fronts couple only through
the total-volume constraint, which does not alter the orthogonal
\(n\ge1\) modes.

**Therefore, for an ideal uniformly extruded channel with neutral
90° plates, the predicted planar local supercritical/subcritical
onset is not pre-empted *at its critical state* by an infinitesimal
vertical instability.** The continuous-versus-subcritical distinction
survives into this rigorously specified 3D capillary problem.

**Limitations:** This is a *local linear stability* statement
and a normal-form implication, not a proof of global stability,
selection of finite-jump landing states, or realistic experimental
plate wetting. It also assumes ideal mobile Young-angle contact lines,
homogeneous sidewall chemistry, no gravity and a regular non-self-
intersecting interface. Finite plate angles other than 90° invalidate
the cylindrical equilibrium and require solving a truly 3D
Young–Laplace boundary-value problem.


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

**Professor assessment:** 5/10 for theoretical readiness; 4/10 for journal submission.
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
