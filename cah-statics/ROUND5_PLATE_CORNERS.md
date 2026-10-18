# Round 5 — Plate-angle perturbation, corner wetting, and a solvability identity

**17–18 October 2026. Model-only computational/theoretical research.**
This round does **not** solve the nonlinear 3D Young–Laplace problem for
arbitrary contact angles, nor prove that the isospectral pair retains
opposite bifurcation classes away from neutral plates.

## A physically important check against Poulose et al. (2023)

Poulose et al., [IEEE Magnetics Letters 14, 3500105](https://doi.org/10.1109/LMAG.2023.3262976),
state that the apparent in-plane angles are about 30°, while the contact
angle at the *upper/lower PMMA plates is close to 90°*. The actual
intrinsic Young angle of the laser-cut sidewall and exact deviation
of the plates from 90° are not independently calibrated in our study.

Therefore the neutral-plate *theorem* of Round 4 is potentially relevant
to the published channel, but the 36.5° angle and exact sinusoidal/
isospectral wall shape remain illustrative parameters, not observations.

## The exact corner gate (existing Concus–Finn theory)

At a geometrically sharp right-angle corner between sidewall and plate,
a smooth interface normal satisfying both Young boundary angles requires

\[
\cos^2\theta_{\mathrm{side}}+\cos^2\theta_{\mathrm{plate}}\le1.
\]

For acute Young angles this is equivalent to

\[
\boxed{\theta_{\mathrm{side}}+\theta_{\mathrm{plate}}\ge90^\circ.}
\]

If the sum is *less* than 90°, the classical **Concus–Finn**
criterion predicts a corner rivulet. This is established prior physics,
**not our discovery** (see the sources below). With the illustrative
sidewall angle 36.5°, plate angles below 53.5° enter the corner-film
regime. With plates close to 90° the threshold is safely avoided in
the ideal sharp-corner model. Exceeding the threshold is **not proof**
that no other film/pinning mechanism exists; rounded corners and
dynamic wetting must also be checked.

In particular **do not** infer from the illustrative 36.5° angle
that the experimental three-layer PMMA chip necessarily has
corner rivulets: its authors explicitly describe plate angles close
to 90°.

## A controlled 3D perturbation, derived and calculated

For each of the two isospectral planar critical cylindrical menisci,
parameterize by planar arclength \(q\in[-L,L]\) and height \(z\in[0,d]\).
At the curvature extremum, the exact normal Jacobi operator is

\[
\mathcal Ju=-u_{qq}-u_{zz}-\kappa_0^2u=\lambda.
\]

The sidewall contact lines are free to move under *linearized*
Young conditions, which supply Robin boundary conditions:

\[
u_q(\pm L,z)=\mp\beta u(\pm L,z), \quad
\beta=\kappa_0\tan\phi,
\]

where \(\phi\) is the circular arc's edge angle measured from its axis.
For a small deviation of the top/bottom intrinsic Young angle from
90°, set \(\epsilon=\cos\theta_{\rm plate}\).
The linearized plate conditions are

\[
u_z(q,0)=-\epsilon,\qquad u_z(q,d)=+\epsilon.
\]

At the planar critical point there is a positive Jacobi zero mode
\(\psi_0(q)=\cos(\kappa_0q)\).
Multiplying the Jacobi equation by \(\psi_0\), integrating, and using
the homogeneous Robin conditions gives the *Fredholm solvability
condition*

\[
\lambda\,d\int_{-L}^L\psi_0\,dq
=-2\epsilon\int_{-L}^L\psi_0\,dq.
\]

The physical curvature perturbation is \(\delta\kappa=-\lambda\)
under our right-meniscus orientation. Consequently

\[
\boxed{\delta\kappa=\frac{2\cos\theta_{\rm plate}}{d}
          +O(\cos^2\theta_{\rm plate})}.
\]

This additive plate-curvature correction is well-known in
rectangular-channel capillarity. Our contribution is a *controlled
critical-state solvability calculation* showing when it follows from
the full capillary Jacobi problem with moving sidewall contact lines,
rather than simply inserting it as a heuristic. **It is the same for
both of our isospectral geometries at first order**: contact-angle
sensitivity of the capillary pressure is not itself a signature of the
Fourier-phase reversal of the quartic stability coefficient.

We solved the *two-dimensional* Jacobi PDE in \(q,z\) by independent
finite differences, simultaneously solving for the pressure multiplier
and imposing \(\int u\,dq\,dz=0\). The last condition chooses a
**zero-mean gauge** along the neutral contact-line translation mode:
reported absolute shape displacements are not unique physical
predictions until the actual volume/driving continuation selects
the translation amplitude.

## Mesh-refined numerical result

Illustrative values are \(w_0=1\) mm, nominal amplitude 0.05 mm,
period 0.3 mm, depth \(d=0.38\) mm, assumed sidewall Young angle
36.5°, and plate Young angle **89°**, just 1° from the analytically
solvable neutral condition.

| Quantity | Minus phase | Plus phase |
|---|---:|---:|
| Neutral critical curvature \(1/R\) (mm⁻¹) | 2.014768 | 1.995566 |
| First-order plate curvature correction (mm⁻¹) | +0.091855 | +0.091855 |
| Relative curvature change | +4.56% | +4.60% |
| Gauge-fixed maximum normal correction at 129×65 mesh | 1.1160 µm | 1.1158 µm |

Three meshes (33×17, 65×33, 129×65) converge in the shape
displacements. The pressure sensitivity agrees with exact
\(2/d=5.2631578947\,\mathrm{mm}^{-1}\) per unit \(\cos\theta_{\rm plate}\).
The first-order curvature contrast from just *one degree* departure
is ~4.6%, larger than the ~1% neutral in-plane pressure contrast
between the two phase-matched examples.

The plotted linear formula beyond ~85° is explicitly **an
extrapolation**, not a validated nonlinear 3D prediction. We do
**not** know the nonneutral quartic coefficient or actual 3D jump
outcomes. The large pressure correction at 1° also motivates testing
the higher-order terms rather than trusting extrapolation.

## Implications for a publishable falsification experiment

1. Measure *separately* the laser-cut sidewalls' intrinsic wetting and
   the upper/lower plates' wetting, with uncertainty. A 30° apparent
   plan-view angle is not necessarily a Young angle.
2. Obtain real wall profiles, including harmonic phases and actual
   corner rounding. Verify 2D YL and finite-depth geometry first.
3. Measure both contact-line positions versus independently estimated
   **liquid volume**, and classify local continuous versus
   discontinuous onset while accounting for pressure/fold competition.
4. Compare against the **published** 2018 and 2024 capillary-snap
   theories and a baseline that only uses roughness power spectra.
5. Numerically continue the *fully nonlinear* 3D free-contact-line
   Young–Laplace surfaces away from 90°; calculate the constrained
   quartic sign and linear modes without assuming the planar values.
   Our new first-order calculation supplies a regression test for that
   more difficult solver.

## Professor's publication assessment

**Journal-submission readiness: 4/10 (unchanged)**. A useful
mathematical physical-domain gate and a reproducible near-neutral 3D
perturbation result have been established. But Concus–Finn is classical,
the first-order plate-curvature term is familiar, and the planar
pitchfork/snap literature is substantial. We still lack a proof that
the *3D nonlinear* phase-sensitive quartic sign survives at realistic
plate wetting, and we have **no blinded quantitative experimental
comparison**. Only after these tests would *Langmuir* or
*Physics of Fluids* become a serious submission target.

## Provenance and replication

- [Poulose et al. 2023, CC BY published article](https://doi.org/10.1109/LMAG.2023.3262976)
- [Bisquert et al. 2024, competing fold/pitchfork snapping](https://doi.org/10.1021/acs.langmuir.4c02908)
- [Spontaneous rise in open rectangular channels 2018](https://pubmed.ncbi.nlm.nih.gov/29793169/)
- [Concus–Finn corner-rivulet criterion review](https://www.nature.com/articles/s41526-021-00180-6)

From the repository root:

```bash
pip install -r cah-statics/requirements.txt
python cah-statics/plate_corner_response.py
```

The script independently generates CSV/JSON diagnostics and PNG/SVG
figures under `cah-statics/outputs/`. No experimental data or videos
were copied, and GitHub Actions is pull-request/manual only.
