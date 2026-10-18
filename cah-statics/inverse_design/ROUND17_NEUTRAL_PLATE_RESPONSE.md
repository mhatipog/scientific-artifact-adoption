# Round 17 — Exact neutral-plate susceptibility and a capillary inverse-identifiability counterexample

**18 October 2026.** Research status: a new exact *model calculation* and a constructive local capillary theorem, with independently convergent FEM and a deliberately retained failed nonlinear mesh gate. **Not an observation, not a validated nonneutral three-dimensional quartic, and not a claim of first-priority novelty.** Existing draft theoretical PR #20 remains unmerged; original Yan optics are isolated in draft PR #14.

## 1. Why this round matters

The previous nominal 170-nm 5/7 Fourier-phase pair has opposite reduced two-front quartics at a selected `32°` sidewall Young angle, but the nonlinear finite-depth discretization had never passed an unbiased convergence test. We subjected a separate Young–Laplace implementation (CasADi/Ipopt) to an exact 90°-plate pressure benchmark. The result is negative: minute KKT residuals coexist with strongly nonconvergent absolute pressure. Numerical stationarity **is not** continuum accuracy.

Instead, an exact continuum neutral-plate Jacobi response can be derived and independently verified, providing a more stringent benchmark and an unexpected separation between linear-response observability and nonlinear two-front instability class.

## 2. Exact fixed-volume plate-wetting susceptibility

Consider one exact cylindrical Young–Laplace meniscus at `theta_plate=90°`, parametrized by arclength `q in [-L,L]` across its circular section, between horizontal plates `z in [0,d]`. The sidewall contact line is mobile at a prescribed intrinsic Young angle `theta_side`. Write the normal displacement `u(q,z)`, curvature `kappa>0`, and support curvature `kappa_wall`. The physical free-contact-line Robin coefficient is

```text
beta = kappa*cot(theta_side) + kappa_wall/sin(theta_side).
```

At the selected pressure-stationary capillary position, `beta = kappa*tan(kappa*L)`, with `0<kappa*L<pi/2`; the nodeless `cos(kappa*q)` is an unconstrained Jacobi zero mode, but **not** a fixed-volume admissible zero mode because its mean is nonzero. The other fixed-volume planar and depth modes are locally positive in the selected neutral-cylinder model.

Let `epsilon=cos(theta_plate)`. The normalized fixed-volume minimized energy (divide physical energy by surface tension `gamma`) admits, within its smooth constrained graph branch,

\[
 E(\epsilon)=E(0)-\frac{2V}{d}\epsilon+Q\epsilon^2+O(\epsilon^3),
 \qquad \frac{dp_{\rm signed}}{d\epsilon}\bigg|_0=-\frac{2}{d}
\]

(the sign of pressure follows the signed energy-volume convention used in the nonlinear solver; its *magnitude* increases by `2 epsilon/d`). The first coefficient is the plate wetting work, not a fitted hydraulic law.

The Jacobi weak equation at first order is

\[
 J=-\partial_q^2-\partial_z^2-\kappa^2,\quad
 J u=\Lambda,\quad
 u_q(\pm L,z)=\mp\beta u(\pm L,z),\quad
 u_z(q,0)=-1,\quad u_z(q,d)=+1,\quad
 \int u\,dq\,dz=0.
\]

Fredholm compatibility with the nodeless zero mode gives `Lambda=-2/d`. Define even normalized Robin eigenfunctions `phi_j(q)=cos(k_j q)/sqrt(N_j)`, with roots `t_j=k_j L` in `(j pi,(j+1/2)pi)` satisfying `t_j tan(t_j)=beta L`, normalization `N_j=L+sin(2t_j)/(2k_j)`, and spatial integrals `G_j=2 sin(t_j)/(k_j sqrt(N_j))`. Set `lambda_j=k_j^2-kappa^2`; `lambda_0=0` and `lambda_j>0` for `j>=1`.

Separation of variables with the volume projection yields the **exact** continuum series

\[
\boxed{
 Q=-\frac{dG_0^2}{6}
 -\sum_{j=1}^{\infty}G_j^2\left[
  \frac{\coth(\sqrt{\lambda_j}d/2)}{\sqrt{\lambda_j}}
  -\frac{2}{d\lambda_j}
 \right] <0.
}
\]

The term inside brackets is `(coth x-1/x)/sqrt(lambda)`, strictly positive. Hence the response is **strictly concave** with respect to plate wetting in this constrained local equilibrium class; this is a model thermodynamic susceptibility, *not* proof of a physical snap.

### A genuine mathematical error bound

For `j>=1`, the root equation implies `|sin(k_j L)|<=beta/k_j`, while `N_j>=L`, `k_j>=j*pi/L`, and `sqrt(lambda_j)>=sqrt(3)*k_j/2`. Also `coth x-1/x<=1` for `x>0`. Therefore truncating after `N` even eigenmodes (indices `0,...,N-1`, `N>=2`) omits a **negative** tail satisfying

\[
0\le Q_N-Q_\infty
\le \frac{2\beta^2L^4}{\sqrt3\,\pi^5(N-1)^4}.
\]

This is an analytical bound on *spectral truncation only*: it is not floating-point interval arithmetic and does not bound uncertainty in measured or assumed physical inputs. We use `N=320`, confirm against `N=512`, and retain the explicit tail bounds.

## 3. Numerical model values and independent finite-element convergence

Nominal geometry: full width `1 mm`, fundamental amplitude `0.05 mm`, period `0.3 mm`, high harmonics `±0.00017[cos(5kx)-cos(7kx)] mm`, depth `d=0.38 mm`, assumed sidewall Young angle `32°`. These are **not fabricated or measured inputs**.

| Quantity | Minus phase | Plus phase |
|---|---:|---:|
| Selected planar critical coordinate `s` (mm) | 0.703953763 | 0.708869057 |
| Circular curvature `kappa` (mm⁻¹) | 2.038783003 | 2.031675054 |
| `Q_320` (mm²) | **−0.087315948234** | **−0.085192444372** |
| Absolute analytical tail bound (mm²) | **1.85e−11** | **9.34e−12** |
| Contrast `Q_plus−Q_minus` (mm²) | — | **+0.002123503862** |

The phase contrast in `Q` is ≈2.46% relative to the mean magnitude (or 2.49% relative to the plus-phase magnitude). This does *not* give a pressure-extremum or snap threshold and cannot determine the fourth-order onset class.

An independent conforming tensor-product P1 FEM assembles `K = (K_q-kappa² M_q)⊗M_z+M_q⊗K_z` plus physical endpoint Robin terms, then solves a **bordered volume-constrained** linear system. The depth mesh spans `[0,d]`—a factor-of-two depth implementation error was discovered and corrected before freezing the tests. The FEM does not call the analytical series.

| Axial × depth nodes | Minus FEM Q (mm²) | Plus FEM Q (mm²) | Approximate error contraction |
|---|---:|---:|---:|
| 17 × 9 | −0.085791106 | −0.083750449 | — |
| 33 × 17 | −0.086931113 | −0.084829600 | ~4 |
| 65 × 33 | −0.087219362 | −0.085101499 | ~4 |
| 129 × 65 | −0.087291766 | −0.085169687 | ~4 |

The largest final absolute error is `2.42e−5 mm²` and decreases about fourfold with mesh spacing halved; constraint and relative stationarity residuals are below `1e−11`. Straight-channel nonlinear CasADi 3D results independently approach the analytic straight-channel value `Q=−0.0753987025 mm²`, albeit at finite `epsilon=cos 89°`.

See [finite-element data](derived/plate_susceptibility_fem.csv), [deterministic SVG](figures/plate_susceptibility_convergence.svg), [series/FEM reproducer](plate_susceptibility.py), [SVG generator](plot_plate_susceptibility.py), and [11 offline tests](tests/test_plate_susceptibility.py).

## 4. A sharper local mathematical obstruction: full linear invisibility of the quartic sign

The exact susceptibility exposes a stronger distinction between **linearized capillary observation** and a **nonlinear symmetry-breaking normal form**.

Take an even channel width `w(x)` and a pressure-stationary meniscus at `x=s>0`. Let `alpha=theta_side+atan(w'(s)/2)`, `z=w'(s)/2`, and `H'=dH/ds>0`. Introduce a mirror-symmetric smooth perturbation around **both** `+s` and `−s` so that the jets of the width through third order vanish at `s`, but `delta w''''(s)=A`. One construction uses an even compact bump `A u^4 chi(u/rho)/24`, with `chi(0)=1` and all derivatives vanishing at its finite support edge. Reflected, nonoverlapping copies preserve channel reflection symmetry; they may also be periodized when supports do not overlap.

The exact circular pressure satisfies

\[
 \frac{\partial p'''(s)}{\partial w''''(s)}
  =-\frac{\sin\alpha}{w(s)(1+z^2)}.
\]

Since `H'`, `H''` and `p''` are unaffected by the fourth jet, the reduced two-front fixed-volume quartic responds **exactly linearly**:

\[
\boxed{\Delta B_4=
\frac{H'\sin\alpha}{12w(1+z^2)}\, A.}
\]

This changes the sign of `B4` at fixed selected critical coordinate and while keeping `w,w',w'',w'''`, pressure/volume linearization, exact neutral-cylinder geometry, *every neutral-plate Jacobi eigenvalue*, and the entire just-derived susceptibility `Q` unchanged. Indeed the first-order plate-forcing solution is identical. The conclusion is exact **within the reduced local two-front model and the neutral-cylinder linearization**.

For a concrete algebraic counterexample, setting `A=-2B4/(dB4/dw'''')` sends `B4 → −B4`:

| Original local phase | Original `B4` | Imposed `delta w''''(s)` (mm⁻³) | New `B4` |
|---|---:|---:|---:|
| Minus | −2413.487 | +119627.779 | +2413.487 |
| Plus | +3605.931 | −175425.905 | −3605.931 |

This is a constructive non-identifiability theorem: no amount of *exact neutral-cylinder linear spectrum or first-order plate-angle response* by itself distinguishes those two quartic classes. Generic bifurcation theory already warns that linearization cannot determine nonlinear coefficients; the contribution here is the explicit Young–Laplace geometric-jet realization and exact sensitivity coefficient, **not a claim to have invented nonlinear bifurcation non-identifiability**.

**Crucial feasibility caveat:** with an illustrative 4-µm bump radius, the maximum width change is only `0.117–0.172 nm`, physically below plausible wall metrology/fabrication control. Letting the radius shrink makes the `C³` disturbance arbitrarily small but concentrates large fourth derivatives. This is a *mathematical counterexample*, **not a robust inverse-designed experimental pair**; it destroys Fourier-power matching of the original 5/7 design, and no global competing-fold exclusion or nonneutral 3D bifurcation continuation follows.

See [constructive code](jet_invisibility.py) and [seven offline tests](tests/test_jet_invisibility.py).

## 5. Nonlinear 3D falsification gate: failed

We independently coded the graph-class full-area surface `r(t,z)=(X(t,z),t*w(X(t,z))/2,z)` with top/bottom plate wetting and sidewall Young wetting, constrained its full liquid volume, and solved using CasADi/Ipopt with tensor P1 geometry and 2×2 Gaussian quadrature. There are no stored video files or Github Actions.

At neutral 90° plates, the pressure **must** equal that of the exact cylindrical solution: `p_minus=−2.03878300318 mm⁻¹`; `p_plus=−2.03167505406 mm⁻¹` with the signed variational convention.

| Axial × depth nodes | Minus P1 pressure | Plus P1 pressure |
|---|---:|---:|
| 17 × 5 | −2.174637 | −2.109304 |
| 33 × 7 | −1.961232 | −2.052614 |
| 49 × 9 | −2.060569 | −2.024494 |
| 65 × 9 | −2.013793 | −2.122198 |
| 97 × 13 | −2.088317 | −1.993799 |

All solver stationarity residuals are tiny, mostly `~1e−15`; nevertheless this pressure is not convergent to its exact value. Likewise a *finite* 90°→89° energy finite difference produces wandering apparent `Q89` values (~−0.064 to −0.080 mm²), not the proven infinitesimal neutral susceptibility. The much more stable pressure increment, approximately `−0.092 mm⁻¹`, reproduces the **already-known** universal leading `−2 cos 89°/d` response but does **not** rescue the unreliable absolute pressure.

Possible causes include nonconvex local branches, insufficient wall resolution and P1 geometric quadrature; **none has been definitively isolated**. Do not extrapolate, Richardson-fit or publish this failing sequence as continuum nonneutral symmetry breaking.

See [optional CasADi pilot reproducer](nonlinear3d_casadi_audit.py) and [rejected ten-row diagnostic](derived/nonlinear3d_failed_mesh_gate.csv). **No nonlinear B4, onset class or global stability assertion is upgraded.**

## 6. Related prior art and novelty restraint

- Hatipogullari et al. (2019), *Physical Review Fluids* 4, 044008, establishes the volume-controlled contact-angle-hysteresis setting (DOI: [10.1103/PhysRevFluids.4.044008](https://doi.org/10.1103/PhysRevFluids.4.044008)).
- Bostwick & Steen (2015), *Annual Review of Fluid Mechanics*, [Stability of Constrained Capillary Surfaces](https://doi.org/10.1146/annurev-fluid-010814-013626), develops established Jacobi/volume/contact-line stability theory. Our eigenfunction Green-response calculation is an application and explicit bound, **not a new Jacobi operator theory**.
- Wells et al. (2018), [*Nature Communications*](https://doi.org/10.1038/s41467-018-03840-6), and Bisquert et al. (2024), [*Langmuir*](https://doi.org/10.1021/acs.langmuir.4c02908), are central prior art on capillary symmetry breaking and fold–pitchfork competition. The Round-17 jet counterexample should be compared to their nonlinear stability frameworks; priority over the broader inverse-spectrum literature has **not** been established.
- Yan et al. (2025) [DOI](https://doi.org/10.1016/j.jcis.2025.137884), a single-front throat–pore experiment, is **not** evidence for this two-front quartic. The original-image physical-depth calibration and shared acquisition-clock questions in empirical PR #14 remain open.

## 7. Reproduction and explicit decision

From `cah-statics/inverse_design`, using NumPy/SciPy for the *mathematical* package:

```sh
python plate_susceptibility.py --csv
python jet_invisibility.py
python plot_plate_susceptibility.py
python -m unittest discover -s tests -p 'test_plate_susceptibility.py' -v
python -m unittest discover -s tests -p 'test_jet_invisibility.py' -v
```

Separately, only when CasADi/Ipopt is installed, the **rejected** nonlinear comparison can be repeated without writing media or dispatching workflows:

```sh
python nonlinear3d_casadi_audit.py --audit
```

Total new passing network-free tests: **18**, in addition to the prior branch suites. The analytical spectral tail has a proved `O(N^-4)` enclosure. No CI workflows were added/triggered manually, no experimental media copied, and the PR should stay **draft and unmerged** until an independent nonneutral 3D solver passes a neutral pressure and energy susceptibility benchmark.

## Professor's assessment

- **Mathematical rigor: 8/10.** Exact Jacobi solution, explicit `O(N^-4)` truncation bound, explicit local fourth-jet theorem and conforming FEM checks; smooth-branch existence and global 3D branches not proven.
- **Novelty relative to literature: 5/10.** The explicit geometric realization is interesting, but the underlying distinction between linear and nonlinear normal forms is established; no robust or physical "first" claim.
- **Experimental validation: 4/10 for programme, 0/10 for this designed pair.** Earlier original Yan data are real, but independent thickness and clock calibration remain unresolved.
- **Potential impact: 7/10.** A falsifiable warning against extrapolating spectral/linear capillary observations into nonlinear hysteresis classifications, plus a certified solver benchmark.
- **Journal-submission readiness: 5/10, unchanged.** More exact mathematics does not supply the missing independently measured prediction.

**Strongest reviewer objection:** An elegant but specialized linear Green-response identity plus a generic high-derivative normal-form counterexample are not yet a novel, experimentally useful law; moreover the only tested corrugated nonlinear 3D solver fails its exact neutral benchmark.

**Possible eventual journal:** *Physical Review Fluids* if a measurable two-front prediction with rigorously converged 3D continuation follows; *Langmuir* if the strongest result becomes experimentally realized wetting-patterned pinning design.

**Single most decisive next calculation:** Implement a geometry-fitted, higher-order free-contact-line nonlinear solver that reproduces **both** the exact 90° absolute pressure and the infinitesimal `Q` above on at least three refinement levels *before* computing the 89° and 88° constrained nonlinear `B4`. Separately, an experimentally calibrated pore-depth map—not another optical intensity threshold—is essential to advance the Yan volume budget.
