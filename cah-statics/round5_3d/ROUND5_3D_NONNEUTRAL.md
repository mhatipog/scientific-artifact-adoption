# Round 5 — Non-neutral finite-depth wetting and an isospectral capillary bifurcation

**17 October 2026. Independently simulated theoretical result, not an experimental comparison.**

## Main result (converged within tested Galerkin families)

For two idealized corrugated channels with precisely the same mean width,
minimum and maximum width, period, RMS height/slope and Fourier power spectrum,
but differing in the phase of a ±0.5 μm harmonic correction, numerical minimization
of a three-dimensional capillary surface confirms *opposite local bifurcation classes*
when top/bottom plate Young angles are **90°, 85° or 80°**.

* Profile A, minus phase: **subcritical** pitchfork (positive third derivative
  of single-front capillary pressure versus displaced volume at the pressure
  extremum).
* Profile B, plus phase: **supercritical** pitchfork (negative third derivative).

The sidewall Young angle is **assumed** 36.5° in both channels. Channel depth is
0.38 mm, nominal width 1.00 mm, period 0.30 mm and peak-to-peak width
variation 0.10 mm. **None of these idealized harmonic microchannels was fabricated
or independently imaged.**

The pressure and bifurcation signs are computed from volume-constrained,
three-dimensional surface-energy minimization, *not* by adding a depth-pressure
correction to the two-dimensional result. We allow the contact lines to move
independently along the top/bottom plates and vertically along the sidewalls,
subject to reflection symmetry in the transverse and vertical directions.

| Sidewall phase | Plate angle | Critical single-front volume (µL) | Critical p=ΔP/γ (mm⁻¹) | Cubic pressure derivative p''' (mm⁻¹/µL³) | Local onset |
|---|---:|---:|---:|---:|---|
| A (minus) | 90° | 0.15514 | 2.01477 | +1623 | Subcritical |
| A (minus) | 85° | 0.15519 | 2.47113 | +1542 | Subcritical |
| A (minus) | 80° | 0.15539 | 2.91935 | +1304 | Subcritical |
| B (plus) | 90° | 0.15800 | 1.99557 | −1761 | Supercritical |
| B (plus) | 85° | 0.15792 | 2.45241 | −1729 | Supercritical |
| B (plus) | 80° | 0.15763 | 2.90200 | −1548 | Supercritical |

These are degree-six local polynomial fits in displaced volume to 19 model
samples spanning ~0.026–0.030 µL, evaluated at the fitted pressure extremum.
The signs agree with degree-four and degree-five fits in this 80°–90° window,
and the 80° values agree in sign under higher surface approximation orders.
The magnitude of p''' is less robust than its sign for the minus phase.

## Surface variational problem

The meniscus is represented as a graph in three-dimensional physical space,
\(\mathbf r=(X(\eta,z), \eta w(X_e(z))/2,z)\), with freely varying
sidewall contact line \(X_e(z)=X(1,z)\). Its free energy is

\[
\mathcal E/\gamma=A_{lg}-\cos\theta_s A_{w,s}
 -\cos\theta_p(A_{w,t}+A_{w,b}),
\]

subject to a fixed displaced liquid volume \(Q\). The volume is computed from
the *actual wall width function* using a cross-sectional Green identity,
rather than assuming the meniscus is cylindrical. Surface-area metric and
four-wall wetting areas are evaluated by Gauss quadrature; derivatives are
obtained by automatic differentiation (JAX). Tensor-product even Legendre
basis modes parameterize the deformation around the **exact planar circular
meniscus**, avoiding spectral errors associated with approximating a
near-tangent cap by an unconstrained low-order polynomial.

The pressure is the negative Lagrange multiplier
\(p(Q)=-d\mathcal E/dQ\). A two-front plug controlled by its total volume
can make one front gain \(\delta Q\) while the other loses \(\delta Q\).
By reflection symmetry, the energy at fixed total volume expands as

\[
\mathcal E(Q+\delta Q)+\mathcal E(Q-\delta Q)
=2\mathcal E(Q)-p'(Q)(\delta Q)^2
-\frac{p'''(Q)}{12}(\delta Q)^4+O((\delta Q)^6).
\]

Thus the symmetry-breaking critical point obeys \(p'(Q_c)=0\). If the other
modes remain stable, the quartic normal-form coefficient has the sign of
\(-p'''(Q_c)\): positive (supercritical) or negative (subcritical).
This clean identity is general for two identical, independently equilibrated
fronts and is *not* presented as a discovery of pitchfork bifurcations.

## Independent numerical controls

1. Neutral (90°) plates: the solver reproduces the exact two-dimensional
   circular-arc Young–Laplace curvature on *both* patterned channels to
   approximately 10⁻¹¹ mm⁻¹. The exact neutral shape is represented as the
   base field; it is **not fitted** in the polynomial basis.
2. Straight rectangular channel: the pressure equals the independent
   force balance \(p=2\cos\theta_s/w+2\cos\theta_p/d\) to roundoff for
   plate angles 90°, 85° and 80°. **That pressure relation is a special
   straight-channel benchmark, not an exact corrugated-channel law.**
3. Surface Galerkin refinement: at 80° the critical-candidate pressure and
   signs are insensitive to moving from 3×2 to 5×4 deformation orders;
   individual pressure differences are ~10⁻⁴ mm⁻¹. The convergence matrices
   and raw traces are preserved as CSV.
4. A separate *parity-complete* constrained Hessian audit includes odd and
   even perturbations across both width and depth (49 and 81 polynomial
   modes), allowing independent sidewall contact-line motion. At the
   representative near-critical states for 80°, 85°, 90°, **all sampled
   constrained eigenvalues are positive**. The values are in coefficient
   coordinates, not continuum physical eigenvalues; only their signs are
   used. This is *not* a mesh-independent continuum stability theorem.
5. The original parameter sweep gave an apparent low-angle (60°–70°)
   bifurcation switch. **We explicitly retract the numerical threshold**:
   its location and even the local sign change fail high-order convergence
   tests. Those raw failure cases are included so future researchers do not
   quietly reintroduce the claim.

## What is *not* established

- No **Fuentes-quality quantitative match** to real jump locations,
  contact-line trajectories, Young angles or capillary-pressure traces.
- No prediction of finite jump duration or unique landing state: local
  stability cannot choose a global path without dynamics/history.
- The equilibria are optimized in the \(y\)- and \(z\)-reflection-symmetric
  sector, with an independent finite all-parity Hessian audit. The full
  nonlinear asymmetric branch and exact continuum Jacobi spectrum are
  not traced across every control volume.
- The physical effect may be smaller than uncontrolled fabrication errors;
  the phase-changing perturbation is 0.5 µm. Real-wall topography must be
  measured more accurately than that.
- The number of modes, contact-line regularity at channel corners,
  and the actual experimental plate wetting remain sources of uncertainty.

## Publication assessment: 5/10

A defensible theoretical/numerical contribution is emerging: a controlled
**isospectral counterexample that survives non-neutral three-dimensional
wetting and finite-depth stability checks**. The application is a prediction
that phase-resolved wall metrology matters more than roughness power spectra
alone for local capillary instability.

However, prior work already studies topography-mediated pitchforks and
competing fold/symmetry-breaking bifurcations. Specifically, Wells et al.
(*Nature Communications*, 2018, doi:10.1038/s41467-018-03840-6) and
Bisquert et al. (*Langmuir*, 2024, doi:10.1021/acs.langmuir.4c02908).
A clean application to *Langmuir* or *Physics of Fluids* would require either
(a) a blind experimental comparison with measured sidewall phases and
wettability, or (b) a strong general theorem/robust dimensional phase map
that outperforms the competing models on existing measurements. As of this
round **the research is not ready for journal submission**, but is a valid
internal scientific result suitable for open reproducibility and discussion.

## Reproduce

Python 3.11+ with `numpy`, `scipy`, `jax[cpu]`, `matplotlib`:

```sh
cd cah-statics/round5_3d
python surface_3d.py                # independent benchmarks and sample surfaces
python scan_3d.py                   # 19-point pressure/volume continuation
python full_stability_3d.py         # finite all-parity Hessian eigenvalue signs
python figures_round5.py            # figures and CSV-based evidence
```

The time-consuming threshold-sensitivity scans are separately supplied under
`research_scripts/` and deliberately excluded from automatic PR runs. The
GitHub Actions **manual/PR-only** workflow executes the lean numerical
smoke tests; there are **no scheduled workflows and no claimed empirical fit**.