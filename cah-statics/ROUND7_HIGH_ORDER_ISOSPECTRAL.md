# Round 7 — Matching all correlations through order five still cannot predict local capillary instability

**18 October 2026 — theoretical counterexample. Not experimentally validated.**

## Stronger than the original isospectral example

For a reflection-symmetric channel, let \(u=kx\) with \(k=2\pi/(0.3\,\mathrm{mm})\), and construct

\[
w_\pm(u)=w_0+a\cos u\pm b[\cos(5u)-\cos(7u)],
\]
\[
w_0=1\ \mathrm{mm},\quad a=0.05\ \mathrm{mm},\quad
b=0.000075\ \mathrm{mm}=75\ \mathrm{nm}.
\]

These profiles are **not translations or reflections** of one another.
The construction enforces:

1. Exactly equal mean and RMS width; exactly equal Fourier magnitude
   spectrum and RMS slope.
2. Exactly equal extrema, \(w_{\min}=0.95\) mm and \(w_{\max}=1.05\) mm.
   Indeed \(|\cos5u-\cos7u|=2|\sin6u\sin u|
   \le12\sin^2u\); because \(b<a/24\), this cannot overtake
   the underlying \(1\pm a\) extrema.
3. **All translation-invariant \(n\)-point correlation functions
   match exactly for \(n=2,3,4,5\).** In Fourier terms every frequency
   tuple \(\nu_1+\cdots+\nu_n=0\), with
   \(\nu_j\in\{\pm1,\pm5,\pm7\}\), has an even number of 5/7 modes
   when \(2\le n\le5\), making its amplitude invariant when \(b\to-b\).
   The 1-point function also matches.

Consequently *all first five centered moments* are equal.
Both distributions have skewness **zero**; equal kurtosis follows
analytically. Thus the earlier Round 3 criticism that skewness
could distinguish the two profiles is resolved in this sharper
counterexample. No claimed equality of sixth moments or all
higher-order correlations.

The first non-invariant frequencies appear at sixth order,
e.g. \((1,1,1,1,1,-5)\).
Writing \(q=w-w_0\), exact Laurent-polynomial calculation gives

\[
\boxed{\langle q_+^6\rangle-\langle q_-^6\rangle
=\frac{3ab}{8}(a^4-10a^2b^2+10b^4).}
\]

Numerically this is \(8.78886474654\times10^{-12}\ \mathrm{mm}^6\).
The wall-profile difference never exceeds approximately
**0.29006 µm**, presenting a genuine nano-scale fabrication and
metrology challenge.

## Yet the local CAH bifurcation classes differ

Apply the exact planar circular Young–Laplace two-meniscus
criterion of Round 2,

\[
B_4 = \frac14\left(H''p''-\frac13H'p'''\right)_{s_*},
\quad p'(s_*)=0.
\]

At an **assumed, not measured** intrinsic sidewall Young angle
\(\theta_Y=36^\circ\):

| Quantity | Minus high-harmonic phase | Plus high-harmonic phase |
|---|---:|---:|
| Critical half-length | 700.19156 µm | 699.45579 µm |
| \(B_4\) (model units) | **−1505.931** | **+812.389** |
| Local normal-form class | Subcritical | Supercritical |

Isolating each pressure stationary root before evaluating the
high-precision derivative criterion, the opposite-sign angular
range is approximately **30.30833° to 37.68774°**, a **7.3794°**
wide theoretical interval.

This result is **strictly a planar nonlinear bifurcation
counterexample**. The full constrained 4-variable Hessian, possible
competing pressure extrema elsewhere, nonlocal meniscus-wall
collisions, and non-neutral plate 3D eigenmodes should also be
independently checked for this new 5/7-harmonic pair. Do not call
these profiles "identical in all roughness properties"; their
sixth-order and higher correlations differ. Do not interpret
local supercriticality as global freedom from snaps.

## Why this moves the scientific argument forward

Round 3 proved PSD insufficiency, but its paired second and fourth
harmonics changed skewness. A reviewer might reasonably object
that a single familiar distributional statistic already distinguishes
the channels. **Here even all correlations through fifth order fail.**

This suggests a precise *information-theoretic* limitation on
finite-order descriptions of geometric roughness as inputs to
deterministic local capillary-instability prediction.
Whether a general theorem for arbitrary order \(N\) can be proved
is an open mathematical research task; one matched example is
not such a theorem. In practical engineering, exact geometry may
be indispensable near narrow curvature-controlled bifurcations.

## Real experiment and practical design trade-offs

The stronger *mathematical* example is not necessarily the most
fabrication-feasible. The earlier 2nd/4th-harmonic phase pair differs
by as much as ~2 µm in wall width (at \(b=0.5\) µm), whereas
this 5th/7th-harmonic construction differs by only ~0.29 µm.
Photolithographic fabrication, calibrated optical/AFM/profilometry
measurements and separately measured Young angles would be needed.
Manufacturing roughness and chemical heterogeneity could dominate
the idealized effect. These are *test proposals*, not completed work.

A useful independent design scan of the original 2nd/4th-harmonic
pair illustrates the trade-off between amplitude and angle tolerance:

| 2nd/4th-harmonic correction | Opposite-sign angle width |
|---:|---:|
| ±0.25 µm | 6.12° |
| ±0.40 µm | 6.63° |
| ±0.50 µm | 6.09° |
| ±0.75 µm | 3.70° |
| ±1.00 µm | 1.91° |

Increasing geometric correction is **not monotonic** in
robustness to wetting uncertainty. The 1 µm case has a much narrower
allowed Young-angle window despite a larger physical difference.

## General finite-order correlation-invisibility lemma

The 5/7 example belongs to a larger construction. For any integer
\(N\ge2\), take an **odd integer \(M>2N\)** and construct

\[
w_\pm(x)=w_0+a\cos(kx)
 \pm b[\cos(Mkx)-\cos((M+2)kx)].
\]

**Claim (proved for the Fourier statistics):** every
translation-invariant \(n\)-point correlation function with
\(1\le n\le N\) is identical for \(w_+\) and \(w_-\),
provided the shapes are analyzed as periodic functions of the
same fundamental period.

**Elementary proof.** A Fourier coefficient occurs only at signed
wavenumbers \(\pm1,\pm M,\pm(M+2)\).
A nonzero term in an \(n\)-point correlation requires a
zero-sum frequency tuple. A product changes sign under
\(b\mapsto-b\) only if it contains an odd number \(r\ge1\)
of high-harmonic factors. The signed high-frequency sum then
has magnitude at least \(M-2r\), because the net coefficient
of the odd-frequency center \(M\) is odd and nonzero.
The \(n-r\) fundamental factors can contribute at most
\(n-r\) in magnitude. But \(M>2N\ge n+r\), hence
\(M-2r>n-r\) and cancellation is impossible.
Thus every zero-frequency product contains an even number
of high harmonics. QED.

The exact maximum/minimum widths are also preserved when
\(b<a/[4(M+1)]\): use
\(|\cos(Mu)-\cos((M+2)u)|
=2|\sin((M+1)u)\sin u|
\le2(M+1)\sin^2u\)
and compare the base-cosine deficit from its extrema.

**What this lemma does NOT prove:** It does not, by itself,
establish that the local bifurcation *always* reverses for
every \(N\). That requires a separate nonzero sensitivity
condition for the fourth-order capillary coefficient; high
derivatives strongly amplify high wavenumbers, but competing
critical roots or loss of physical admissibility can intervene.

### A higher-order numerical illustration

For \(N=7\), choose \(M=15>14\) and
\(b=0.046875/M^4\) mm, equivalent to
\(b\approx0.92593\) **nm**. At the unperturbed critical
angle \(\theta_Y=37.14514644^\circ\), direct high-precision
derivatives and single critical-root isolation produce
\(B_{4,-}\approx-932.54\) and \(B_{4,+}\approx+866.50\).
The lemma guarantees exact equality of correlations **through
seventh order**. We have not checked the nonlinear 3D menisci
or global dynamics in this constructed case.

These nanometre perturbations are clearly beyond what we
would responsibly treat as a validated continuum microfluidics
fabrication prediction. The point is a **mathematical
identifiability limitation**, potentially generalizable to a
rigorous theorem with independently established capillary
nondegeneracy assumptions.


## Priority and prior art

The strongest immediate next tasks, in order, are:

1. Independently verify the reduced-model 4-DOF constrained Hessian,
   wall clearance and all nearby branches for this exact 5/7 example.
2. Extract **three-dimensional** \(B_4\) by fixed-volume continuation
   with mesh and solver-order convergence, using a bordered-Hessian
   / Lyapunov–Schmidt method rather than fitting noisy fourth powers.
3. Compare with Wells et al. 2018 and Bisquert et al. 2024;
   their smooth-topography fold/pitchfork physics predates our
   investigations. Verify that this *higher-order-statistical
   degeneracy* is not already in wetting literature.
4. Only then test measured wall profiles (blinded) or fabricate
   controlled isospectral geometries.

[Wells et al. (2018)](https://doi.org/10.1038/s41467-018-03840-6);
[Bisquert et al. (2024)](https://doi.org/10.1021/acs.langmuir.4c02908).

### Independent analytic extraction approach

For a scalar symmetry-breaking coordinate \(c\) and stable orthogonal
degrees of freedom \(y\), write the constrained capillary energy
around the critical symmetric equilibrium as
\[
F(c,y)=F_0 +\tfrac12 y^T H y +\tfrac12 c^2 g^T y
 +\tfrac1{24}F_{cccc}c^4+\cdots.
\]
Minimizing in \(y\) gives the **exact local quartic coefficient**
\[
\boxed{B_4=\tfrac1{24}F_{cccc}-\tfrac18 g^T H^{-1}g},
\quad g=F_{ycc}.
\]
This is the Lyapunov–Schmidt / Schur-complement reduction of a
nondegenerate constrained variational problem, **not a new theorem**.
Apply it after projecting onto the fixed-volume tangent subspace,
isolating the critical zero mode, and checking that other modes
have positive stiffness. It avoids unstable finite-difference
estimates of tiny fourth-order energy differences.
For symmetry-preserving finite-depth continuation,
neutral-plate \(B_4\) should recover the 2D result after consistent
energy/coordinate normalization.

## Professor's assessment — Round 7

Theoretical insight **8/10**;
novelty confidence **5/10** (priority not yet established);
experimental validation **1/10**;
journal submission readiness **4/10**, unchanged.

This is potentially a cleaner publishable mathematical finding,
but not a publishable capillary-physics paper yet. Journal targets
would remain *Langmuir* or *Physics of Fluids* after the 3D and
experimental validation. *Physical Review Fluids* needs a broadly
relevant generalized theorem or strong physical demonstration.

Reproduction:
```bash
python cah-statics/spectrum_phase_sixthorder.py
```
Outputs under `cah-statics/outputs/`. Theory only.
