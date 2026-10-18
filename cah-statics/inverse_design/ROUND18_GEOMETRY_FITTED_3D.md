# Round 18 — A geometry-fitted nonlinear 3D Young–Laplace benchmark finally converges

**18 October 2026. Status: accepted local nonlinear numerical benchmark, conditional model only.** This round removes a concrete *near-neutral, single-meniscus* finite-depth convergence obstacle. It **does not prove a nonneutral two-front quartic sign, global capillary stability, physical snap distance, experimental Young angle, or observed phase contrast.** Draft PR #20 remains unmerged; Yan original-image volume calibration remains a separate unresolved draft PR #14.

## 1. What failed before—and the audited remedy

Round 17 introduced an independent CasADi full-area variational 3D graph solver but **rejected its corrugated neutral-plate pressures**: coarse/refined grids with small KKT residuals missed the known cylindrical pressure, sometimes by 0.05–0.15 mm⁻¹. Simply increasing uniform transverse resolution did not give orderly convergence; unconstrained Ipopt sometimes followed distant branches or even folded parameterizations. Those invalid sequences remain preserved in [the previous report](ROUND17_NEUTRAL_PLATE_RESPONSE.md).

This round independently tested:
1. **Exact cylinder embedding** at every quadrature point, rather than interpolating its high-curvature geometry with piecewise-linear coordinates.
2. Uniform transverse and Chebyshev-clustered meshes; neither reliably solved both phases at moderate cost.
3. A **physical-axial-distance mesh** and a local, volume-constrained stationarity solver instead of global energy minimization.

The successful mesh is derived analytically. For the nominal stationary neutral cylinder, let
\[
R=\frac{w(s)}{2\cos\alpha},\quad
x_c=s+R\sin\alpha,\quad
x_{\mathrm{tip}}=x_c-R,\quad
\alpha=\theta_Y+\arctan[w'(s)/2].
\]
Choose positive-side axial nodes `x_j` uniformly in `[x_tip,s]`, and invert the exact cylinder in the wall-fitted transverse coordinate:
\[
t_j=\frac{2\sqrt{R^2-(x_j-x_c)^2}}{w(x_j)}\in[0,1],
\]
then reflect to `[-1,1]`. Verify strict monotonicity; reject a geometry that violates it.

The full area—not its planar approximation—is integrated for
\[
\mathbf r(t,z)=(X(t,z),\ t w(X(t,z))/2,\ z),
\quad X=X_{\mathrm{cylinder}}(t)+u_h(t,z).
\]
The fixed-volume liquid energy includes the sidewall solid-wetting primitive and the two plate wetting areas:
\[
\mathcal E_\epsilon[u]=
\int\sqrt{X_t^2+[w(X)/2+t w'(X)X_t/2]^2+[w(X)X_z/2]^2}\,dt\,dz
-\cos\theta_Y\sum_{\rm walls}\int S(X_{\rm wall})dz
-\epsilon\,A_{\rm plate,wet},
\]
with `S'(x)=sqrt(1+(w'(x)/2)^2)`, `epsilon=cos(theta_plate)`, and the exact variable-domain liquid volume integral. A local sparse Newton solve enforces `grad(E-pV)=0` **and the prescribed volume** simultaneously.

All retained stationary solutions pass post-solve gates for KKT and volume residuals, small departure from the embedded exact neutral cylinder, and positive `dy/dt=[w+t w' X_t]/2` at every element quadrature point. **Positive discrete constrained Hessian eigenvalues were separately checked** on two meshes. No global search is used or claimed.

**Important engineering lesson:** the uniform transverse coordinate compresses physical Fourier oscillations into very small intervals in the minus phase (its `max |dX_cylinder/dt|≈4.985 mm`, versus `2.398 mm` for plus). Physical-distance grading, rather than merely increasing an arbitrary node count, addresses the true resolution scale. This is a numerical remedy, not a claim of a novel mesh-generation principle.

## 2. Absolute neutral-pressure benchmark (not fitted)

Nominal channel geometry and sidewall angle are the same hypothetical 5/7 phase pair as Rounds 14–17: `d=0.38 mm`, `theta_Y=32°`, width `1+0.05 cos(kx)±0.00017(cos(5kx)-cos(7kx)) mm`. The values **are not independently measured**.

The exact 90°-plate signed pressures are `p_minus=−2.03878300318053 mm⁻¹` and `p_plus=−2.03167505405863 mm⁻¹`.

At 257 planar × 17 depth nodes, the locally solved nonlinear graph gives signed-pressure errors:
- **Minus:** `+9.10207×10⁻⁶ mm⁻¹` (relative `4.46×10⁻⁶`).
- **Plus:** `+3.53914×10⁻⁶ mm⁻¹` (relative `1.74×10⁻⁶`).

At 513 × 3 nodes, the *neutral cylinder's planar pressure alone* has residual ≈`1.27×10⁻⁶ / 7.02×10⁻⁷ mm⁻¹` (minus/plus). This smaller depth mesh is **not** itself validation of a depth-varying nonneutral interface.

Across tested 65–257 × 5–17 meshes, neutral errors decrease in the finer levels; the previously oscillatory uniform grid is **not** retroactively validated.

## 3. Independent finite-depth nonneutral energy benchmark

The exactly derived Round-17 local fixed-volume coefficient is
\[
E(\epsilon)=E(0)-2\epsilon V/d+Q_0\epsilon^2+O(\epsilon^3),
\qquad \epsilon=\cos\theta_{\rm plate}.
\]
The `Q_0` references are independent continuum Jacobi spectral solutions with a proved `O(N^{-4})` analytical tail bound:
`Q_{0,-}=−0.0873159482495` and `Q_{0,+}=−0.0851924443800 mm²`.

We solved the **full nonlinear free-depth** stationary graph at `theta_plate=89.9°` using exactly the same nominal material and geometry inputs, and calculated the finite-angle
\[
Q_{\epsilon,h}=\frac{E_h(\epsilon)-E_h(0)+2\epsilon V_h/d}{\epsilon^2}
\quad\text{(not a direct measurement of the exact derivative at zero)}.
\]

| Physical-axial × depth nodes | Minus nonlinear Q (mm²) | Plus nonlinear Q (mm²) |
|---|---:|---:|
| 65 × 5 | −0.08124305 | −0.07922200 |
| 97 × 7 | −0.08401889 | −0.08253268 |
| 129 × 9 | −0.08540989 | −0.08368900 |
| 193 × 13 | −0.08645714 | −0.08452054 |
| 257 × 17 | **−0.08682926** | **−0.08481289** |
| Independent neutral continuum Q | **−0.08731595** | **−0.08519244** |

Finest raw relative differences are **0.557% (minus)** and **0.446% (plus)**. All reported fixed-volume KKT residuals are near double precision, with physical graph displacements roughly 0.6–1.2 µm at 89.9°; transverse Jacobians remain positive.

### Independent axial/depth refinement and a posteriori Richardson check

We separately changed axial and depth discretization at the following *four paired meshes*, holding the wetting assumptions, physical geometry and small angle fixed:

| Phase | Q(129,9) | Q(129,17) | Q(257,9) | Q(257,17) |
|---|---:|---:|---:|---:|
| − | −0.0854098912 | −0.0864728298 | −0.0857674739 | −0.0868292621 |
| + | −0.0836890005 | −0.0847159477 | −0.0837861199 | −0.0848128869 |

As a check—not a proof or experimental calibration—we remove assumed leading `O(h_x²+h_z²)` terms by
\[
Q_{\mathrm{R}}=
\frac{16Q_{257,17}-4Q_{129,17}-4Q_{257,9}+Q_{129,9}}{9}.
\]
The predicted continuum comparisons are:

| Phase | Double-extrapolated nonlinear Q (mm²) | Independently exact Q (mm²) | Relative difference |
|---|---:|---:|---:|
| − | −0.0873018744 | −0.08731594825 | **0.0161%** |
| + | −0.0851874356 | −0.08519244438 | **0.00588%** |

The analytical phase contrast is `Q_plus−Q_minus≈0.00212350387 mm²`; the doubled Richardson estimate is `0.00211443882 mm²`, within **0.427%** of that *model* phase contrast.

**Do not mistake these extrapolated accuracies for guaranteed continuum tolerances.** The second-order ansatz was assessed retrospectively against the exact reference, and residual truncation/cross-terms have not been rigorously bounded. Nevertheless, independently changing each grid direction and matching a completely different Jacobi continuum formulation is strong corroborating numerical evidence.

An independent **linear** conforming P1 Jacobi FEM reproducer further isolates depth truncation. At axial 129 nodes, its Q values for minus are `−0.0816853 (z=5)`, `−0.0858967 (z=9)`, `−0.0869586 (z=17) mm²`. The difference between these and the nonlinear values on the 129-node physically fitted axial mesh stays near `+0.000487 mm²` across these three depth grids. For plus, the analogous difference stays near `+0.000132 mm²`. At 257 axial nodes these residuals decrease to `+0.000131` and `+0.000036 mm²`, consistent with a *separable leading axial discretization bias*. This is **numerical evidence of error structure**, not a general error theorem.

### Finite-angle check

On a fixed 129 × 9 mesh, repeating the nonlinear solve at plate angles `89.95°, 89.9°, 89.8°` shifts the apparent Q by only a few `10⁻⁶ mm²` per 0.05–0.1° increment. The signed pressure changes obey the pre-existing universal first-order law `dp/d(epsilon)=−2/d` to the shown precision. No nuisance parameter was fitted.

## 4. Constrained stability: numerical, not continuum

The smallest numerically projected constrained KKT Hessian eigenvalues are strictly positive:

| Phase | 65 × 5, 90° | 65 × 5, 89.9° | 129 × 9, 90° | 129 × 9, 89.9° |
|---|---:|---:|---:|---:|
| − | 0.0012490 | 0.0012705 | 0.0003525 | 0.0003528 |
| + | 0.0026981 | 0.0026999 | 0.0005988 | 0.0005990 |

The values depend on the discrete nodal-coordinate normalization and **should not** be compared as continuum eigenvalues across grids. Their **positive signs at each tested mesh** exclude discrete constrained negative modes *within the represented graph class* for those local states. This still does **not** establish continuum Hessian positivity at nonneutral plates, nor global energy minimization, nor two-meniscus quartic stability.

## 5. What is and is not now established

**Established analytically:** exact 90° cylinder and Round-17 Jacobi susceptibility at local critical menisci; proved series truncation bound.

**Now substantially stronger numerical evidence:** independent 3D nonlinear variational solver, physical-distance meshing, local KKT continuation, exact-neutral-pressure recovery, fixed-volume and orientation checks, positive sampled constrained discrete Hessian, two-axis refinement, and agreement of near-neutral Q with independent linear Jacobi continuum to <1% on the finest raw meshes and <0.02% on exploratory extrapolation.

**Not established:** nonlinear three-dimensional two-front `B4` sign at nonneutral plates; competing global capillary folds, folds near a physical volume trajectory, experimental contact-angle hysteresis or jump distance; measured wall Young angle; microscope-derived liquid thickness; shared camera–pressure logger clock; originality priority for numerical mesh grading. The paired Fourier spectra are identical **only** for the ideal nominal CAD profiles.

Related established prior art includes Hatipogullari et al. (2019), [DOI](https://doi.org/10.1103/PhysRevFluids.4.044008); Bostwick & Steen (2015), [stability review](https://doi.org/10.1146/annurev-fluid-010814-013626); Wells et al. (2018), [DOI](https://doi.org/10.1038/s41467-018-03840-6); Bisquert et al. (2024), [DOI](https://doi.org/10.1021/acs.langmuir.4c02908); and established volume-constrained capillary surface minimization and mesh techniques. **No priority claim** for arc-length mesh refinement or nonlinear KKT continuation is made.

## 6. Reproduction, computational costs and decision

- [Local free-contact-line 3D code](nonlinear3d_local_continuation.py) — optional CasADi, NumPy, SciPy; no JAX.
- [Seven new offline regression tests](tests/test_nonlinear3d_local_continuation.py) — tested successfully, including graph orientation, neutral pressure, nonneutral Q, angle halving and discrete Hessian.
- [10-row diagonal mesh comparison](derived/nonlinear3d_local_diagonal.csv).
- [Independent axial/depth refinement record](derived/nonlinear3d_local_grid_axes.csv) and [2×2 Richardson check](derived/nonlinear3d_local_richardson.csv).
- [Deterministic SVG of raw mesh error](figures/nonlinear3d_local_convergence.svg), with [generator](plot_nonlinear3d_local.py).

Run with an existing local CasADi installation and **no GitHub workflow**:

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_nonlinear3d_local_continuation.py' -v
python cah-statics/inverse_design/nonlinear3d_local_continuation.py --fast --csv
# Full ten-row audit is an optional more expensive CPU run:
python cah-statics/inverse_design/nonlinear3d_local_continuation.py --csv
# Independent two-axis check is also optional and expensive:
python cah-statics/inverse_design/nonlinear3d_local_continuation.py --grid --csv
python cah-statics/inverse_design/plot_nonlinear3d_local.py
```

No experimental videos, numerical fits to measured pressure, workflows or runner dispatches were added. The previous **failed unconstrained** solver and its negative data remain available; they must not be overwritten by this success. The draft branch is **not ready to merge for journal claims** without cross-solver verification of nonneutral two-front bifurcation coefficients.

### Professor's assessment

- **Mathematical rigor 8/10:** Round-17 exact continuum Jacobi series and fourth-jet identity remain the proofs; this round adds stronger but nonrigorous 3D numerical convergence.
- **Novelty 5/10:** This is a meaningful recovery of a difficult benchmark in our geometry, but physical-distance meshing, KKT continuation and Richardson extrapolation are established numerical methods. No independent literature-priority finding.
- **Experimental validation 4/10 for the programme, 0/10 for the phase pair:** No new measured liquid-volume or wettability input.
- **Potential scientific impact 7/10:** Nonneutral 3D numerical predictions can now be checked against an exact independent benchmark, enabling the next critical bifurcation test.
- **Journal-submission readiness 5/10, unchanged:** Progress in solver reliability is not an independently verified new physical prediction.

**Strongest reviewer objection:** all quantitative agreement is with the programme's own analytical model using hypothetical geometry and wetting; a real two-meniscus nonneutral bifurcation type and measured pressure–volume path remain unavailable.

**Potential journal:** *Physical Review Fluids* contingent upon an independent nonneutral 3D two-front normal-form coefficient and a measured phase-engineered channel or equivalent independent experimental pressure–volume validation.

**Single decisive next calculation:** Use this now-benchmarked geometry-fitted KKT solver, plus an independent second implementation, to compute the *two-front constrained quartic* `B4` at 90°, 89.9°, 89° and 88° with full Hessian and all competing folds checked. Do not infer the quartic sign from the validated single-meniscus `Q` alone.

**Remaining experimental blocker:** original Yan pore-filling images still require calibrated incremental water depth and clock provenance; the optical intensity proxy cannot be promoted to physical volume. Treat Yan single-front evidence independently from the two-front design.
