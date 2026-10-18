# Round 6: nonlinear 3D equilibria and first-order invariance

**18 October 2026. Scientific status:** We obtained constrained stationary
menisci from a **fully nonlinear surface-graph energy**, with freely sliding
sidewall and top/bottom contact lines, at prescribed liquid volume.
The calculation does **not** determine the nonlinear three-dimensional
quartic bifurcation coefficient, distinguish full 3D snap trajectories,
or match an experimental video.

## 1. Exact geometric energy, not a guessed 3D curvature

Parameterize the right liquid/gas interface as

\[
{\bf r}(t,z) =
  \left(X(t,z),\tfrac12 t\,w[X(t,z)],z\right),
\quad -1\le t\le1,\quad 0\le z\le d.
\]

Here \(X(t,z)\) is an unknown nonlinear field, including both wall
contact-line curves \(X(\pm1,z)\) and the top/bottom traces.
Let \(B'(x)=w(x)\), \(S'(x)=\sqrt{1+w'(x)^2/4}\).
The exact volume and total (surface-tension-normalized) energy for
a **single right meniscus** in this graph class are

\[
V[X]=\tfrac12\int_0^d\int_{-1}^1B(X)\,dt\,dz,
\]

\[
\begin{split}
E[X;\epsilon]&=\int_0^d\int_{-1}^1
\sqrt{X_t^2+
       [w(X)/2+t\,w'(X)X_t/2]^2+
       [w(X)X_z/2]^2}\,dt\,dz\\
&\quad-\cos\theta_s\int_0^d
   [S(X(-1,z))+S(X(1,z))]\,dz\\
&\quad-\frac{\epsilon}{2}\int_{-1}^1
   [B(X(t,0))+B(X(t,d))]\,dt,
\end{split}
\]

where \(\epsilon=\cos\theta_{\mathrm{plate}}\).
Irrelevant constant reference wetting areas have been omitted.
We minimize \(E\) at fixed \(V\), differentiating the exact
finite-element quadrature with JAX automatic differentiation;
SLSQP handles volume conservation. This is a genuine *nonlinear
3D equilibrium solve within a single-valued graph class*, not
a fully general non-graph meniscus/free-topology method.

The corrugated pair has mean width 1 mm, wavelength 0.30 mm,
fundamental amplitude 0.05 mm, harmonic phase corrections
\(\pm0.0005[\cos(2kx)-\cos(4kx)]\) mm, depth 0.38 mm,
and **assumed** intrinsic sidewall Young angle 36.5°.
These are model inputs, not measured chip parameters.

## 2. Independent exact validation

In a straight uniform rectangular channel, force balance gives

\[
\frac{|\Delta P|}{\gamma}
 =\frac{2\cos\theta_s}{w}
 +\frac{2\cos\theta_{\rm plate}}{d}.
\]

Our nonlinear solver reproduces the exact pressure at plate
angles **90°, 89°, and 87°** on a 25×13 finite-element grid
to within \(3\times10^{-4}\,\mathrm{mm}^{-1}\).
The formula itself is established capillarity, **not a novel result**.

## 3. A sharper result than the pressure shift

At neutral plates (\(\epsilon=0\)), all cylindrical menisci
have identical liquid cross sections at both plates and every
intermediate depth. Hence, for **every** cylindrical trial
configuration at total fixed liquid volume \(V\),

\[
A_{\rm wet,top+bottom}=\frac{2V}{d}.
\]

The plate-wetting perturbation therefore adds exactly
\(-2\epsilon V/d\) at first order, independent of plug centre,
meniscus separation, or Fourier phases. With the envelope theorem
and the usual smoothness/nondegeneracy assumptions after fixing
the bifurcation coordinate \(c\), the *relaxed* energy satisfies

\[
\boxed{E_{\rm eff}(V,c;\epsilon)=E_{\rm eff}(V,c;0)
                 -2\epsilon V/d+O(\epsilon^2).}
\]

It follows that at fixed V the translation-mode quadratic and quartic
coefficients have no linear-in-\(\epsilon\) correction:

\[
\boxed{B_4(\epsilon)=B_4(0)+O(\epsilon^2)}
\]

along a properly continued local critical branch. Thus the distinct
**nonzero** neutral-plate quartic signs from Round 3 persist in
some (not yet numerically bounded) sufficiently small neighborhood
of 90° under branch regularity, but **we have not computed a 3D
\(B_4\) at 89° or any other non-neutral angle**.

This is consistent with the **first-order** pressure-curvature
correction \(2\epsilon/d\): an energy term linear in V changes
pressure but not the instability landscape at fixed V.

The leading nonlinear shape effect is quadratic in \(\epsilon\).
For a locally stable constrained branch, the standard Hessian/envelope
formula is \(-\epsilon^2\langle g_\perp,L_\perp^{-1}g_\perp\rangle/2\),
with projected plate forcing \(g_\perp\). This term can distinguish
different wall phases.

## 4. Numerical results and honest convergence limitations

For \(\theta_{\rm plate}=89^\circ\), the first-order theory gives
a curvature increment of **0.09185477 mm⁻¹**.
Our locally executed nonlinear equilibria give:

| transverse×depth mesh | Phase A (negative correction) | Phase B (positive correction) |
|---|---:|---:|
| 33×9 | 0.094507 mm⁻¹ | 0.091992 mm⁻¹ |
| 65×9 | 0.091758 mm⁻¹ | 0.092026 mm⁻¹ |
| 97×9 | 0.093235 mm⁻¹ | 0.091522 mm⁻¹ |

**Not fully converged:** the absolute corrugated-channel pressures
vary more than these *differences* when refining the transverse
mesh, and even the increment errors are nonmonotone. These results
support a leading-order prediction at the roughly few-percent level;
they do **not** justify precision 3D snap thresholds or quoted
finite-jump lengths.

A more robust *energy* diagnostic uses
\(E_{\rm 3D}(\epsilon)-E_{\rm 3D}(0)+2\epsilon V/d\).
For the 65×9 mesh, dividing by \(\epsilon^2\) gives:

| plate angle | Phase A | Phase B |
|---|---:|---:|
| 89.5° | −0.079843 | −0.075317 |
| 89° | −0.079826 | −0.075359 |
| 88° | −0.079826 | −0.075492 |

The approximately constant values are independent numerical
support for **second-order** shape sensitivity. Values are
normalized model energies, not measured physical energies.
The calculation is constrained to stable *single-front*
equilibria; it does not follow both fronts through a 3D
pitchfork or compute a 3D stability spectrum.

All runs reported constrained optimizer success and volume error
below \(10^{-10}\,\mathrm{mm}^3\). Optimization residuals were
below \(10^{-5}\) in finite-element gradient coordinates. The
straight-channel exact pressure test independently checks
the functional and contact-line variations.

## 5. Novelty and publishability audit

Existing published work already covers patterned-surface snapping,
fold/pitchfork competition, and general capillary variational
principles:

- [Wells et al. 2018, Nature Communications](https://doi.org/10.1038/s41467-018-03840-6).
- [Bisquert et al. 2024, Langmuir](https://doi.org/10.1021/acs.langmuir.4c02908).
- [Langmuir 2023: rectangular channels with varying wetting](https://doi.org/10.1021/acs.langmuir.3c02865).
- [oomph-lib: free-contact-line Young–Laplace FEM reference](https://oomph-lib.github.io/oomph-lib/doc/young_laplace/contact_angle/html/index.html).

**Professor interpretation:** Useful, defensible progress in nonlinear
physical verification. The first-order invariance follows standard
variational perturbation theory and is unlikely to stand alone as
a novel journal contribution. The phase-only roughness
counterexample *may* provide a new design/metrology insight,
but nonlinear **3D quartic coefficients and experimental data are
still missing**. Readiness **4/10**; potential journals
*Langmuir* or *Physics of Fluids* only after those strict gates.

The decisive next work is: adaptively mesh the near-grazing sidewall
contact regions; continue a **two-front** 3D constrained branch;
extract independent 3D Hessian and quartic coefficients, and compare
blindly with a measured channel/video sequence.

## Replication and runner budget

```bash
pip install -r cah-statics/requirements-3d.txt
python cah-statics/nonlinear3d_variational.py
python cah-statics/nonlinear3d_variational.py --research
```

The first run does the exact straight-channel benchmark and 33×9/65×9
corrugated tests. The optional flag also runs the slower 97×9 test.
Results are in `cah-statics/outputs/`; the new workflow is **manual
only**, not automatically scheduled or run on every PR.
