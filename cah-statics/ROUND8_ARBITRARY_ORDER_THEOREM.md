# Round 8 — Arbitrary finite-order correlation invisibility and capillary bifurcation reversal

**18 October 2026 · Scope: rigorous Fourier argument plus a perturbative
local Young–Laplace bifurcation theorem, numerical validations, and an
exact special-case neutral-plate 3D lift.**

**Never extrapolate these examples to fabricated subatomic roughness or
claim a non-neutral 3D pitchfork measurement.**

## 1. Quantified theorem

Fix an integer \(N\ge2\). Choose an odd \(M>2N\) and let

\[
w_{\pm,M}(x)=w_0+a\cos(kx)\ \pm\
     \frac{\beta}{M^4}
     [\cos(Mkx)-\cos((M+2)kx)]
\]

with \(w_0=1\) mm, \(a=0.05\) mm and \(k=2\pi/(0.3\,{\rm mm})\).
We use \(\beta=0.046875\,{\rm mm}\) in the numerical examples.

**Theorem A (exact roughness statistics).** All translation-invariant
\(n\)-point correlation functions, at arbitrary relative offsets, agree
for the \(+\) and \(-\) profiles for every \(1\le n\le N\).
For all sufficiently large \(M\), their mean, RMS, Fourier power
spectra, and extrema \(w_{\min}=0.95\) mm,
\(w_{\max}=1.05\) mm also coincide exactly. They are not related
by a translation or reflection.

*Proof.* Each centered height mode has wavevector in
\(\{\pm1,\pm M,\pm(M+2)\}\). Only mode products whose wavevectors sum
to zero survive spatial averaging. Any product that changes sign
under \(\beta\mapsto-\beta\) contains an odd number \(r\) of the
high-frequency factors. Writing the signed high-frequency sum as
\(\ell M+2m\), \(\ell\) is odd and
\(|\ell M+2m|\ge M-2r\). The remaining \(n-r\) fundamental factors
have total magnitude at most \(n-r\). Since \(M>2N\ge n+r\), they
cannot cancel. All zero-sum terms thus have even \(r\) and identical
coefficients. The extrema follow from
\(|\cos Mu-\cos(M+2)u|
 =2|\sin((M+1)u)\sin u|
 \le2(M+1)\sin^2u\): it suffices that
\(\beta/M^4<a/[4(M+1)]\). QED.

This is **a direct variation of the general phase-retrieval /
homometry phenomenon**, not the first discovery of nonuniqueness of
spectral data. See Grimm and Baake, *Homometric Point Sets and Inverse
Problems* (2008): https://arxiv.org/abs/0808.0094.

## 2. Asymptotic Young–Laplace bifurcation reversal

Use the exact planar two-front capillary model with fixed volume, wall
pressure \(p(x)=2\cos[\theta+\arctan(w'(x)/2)]/w(x)\), and volume
primitive \(H(x)\) defined in the Round 2 report. The local pitchfork
point \(s_*\) satisfies \(p'(s_*)=0\), \(p''(s_*)\ne0\), and the quartic
coefficient is

\[
B_4=\frac14\left(H''p''-\frac13H'p'''\right)_{s_*}.
\]

The sinusoidal model at \(\theta_*=37.1451464361036^\circ\)
has a **nondegenerate** pitchfork with \(B_4=0\) at

\[
s_0=0.698614069163614\ {\rm mm},\quad
p_0''(s_0)=-139.0547168314,\quad
H_0'(s_0)=3.02224362445.
\]

This numerical root is **not an interval-certified proof of exact
parameter identity**; the theorem can instead be phrased for the
nearby exact roots defined implicitly by the analytic equations and
regularity conditions.

**Theorem B (local capillary reversal for every finite-order
statistical limit).** At that simple planar local-degeneracy point,
with \(\beta>0\) fixed and \(M\to\infty\) odd,

\[
\boxed{
B_{4,\pm}(M)=\pm\beta\mathcal C\{
\cos(Mu_0)-\cos[(M+2)u_0]\}
+O(M^{-1}),
}
\]

where \(u_0=ks_0\) and

\[
\mathcal C=
\left[
\frac{H'\sin(\theta+\arctan(w'/2))k^4}
{12w(1+(w'/2)^2)}
\right]_{w=w_0+a\cos(kx),\ x=s_0}.
\]

Our independent high-precision evaluation gives
\(\mathcal C=8793.28959470\) in the prescribed units.
Since \(u_0/\pi=4.65742712776\) is not a multiple of \(1/2\),
\[
\cos(Mu_0)-\cos((M+2)u_0)
  =2\sin((M+1)u_0)\sin u_0
\]
is bounded away from zero on an infinite subsequence of odd \(M\).
Thus for **every** finite \(N\) one can choose
arbitrarily high odd \(M>2N\) such that the two local quartic
coefficients have opposite signs.

*Perturbative proof sketch, with controlled orders.* For
\(b=\pm\beta/M^4\), perturbations of \(w,w',w'',w'''\)
are \(O(M^{-4}),O(M^{-3}),O(M^{-2}),O(M^{-1})\).
The shift in the isolated \(p'\) root is \(O(M^{-2})\),
using \(p_0''\ne0\) and the implicit-function theorem.
Every contribution to \(B_4\) vanishes at least as \(O(M^{-1})\)
except the high-frequency contribution to \(p'''\) through \(w''''\).
Indeed

\[
\frac{\partial p}{\partial w'}
=-\frac{\sin[\theta+\arctan(w'/2)]}
{w(1+(w'/2)^2)},
\]

so the fourth derivative of the small wall disturbance generates
the boxed leading term. The \(O(M^{-2})\) contact displacement
changes its high-frequency phase by only \(O(M^{-1})\).
Higher nonlinear combinations are also lower order because
the second and third derivatives of the disturbance vanish.
The phase factor is nonzero infinitely often, giving opposite signs.
This is a local static normal-form statement, **not an experimental
jump-selection theorem**.

### High-precision sample

| M | |b| (nm) | B4 minus | B4 plus |
|---:|---:|---:|---:|
| 15 | 0.925926 | -932.54 | +866.50 |
| 25 | 0.120000 | +414.92 | -410.33 |
| 51 | 0.006929 | -494.81 | +506.15 |
| 75 | 0.001481 | +29.07 | -33.03 |
| 101 | 0.000450 | +176.02 | -174.20 |
| 151 | 0.000090 | +137.81 | -140.35 |
| 201 | 0.0000287 | -420.28 | +416.55 |

Both signs reverse *relative to each other*, not necessarily with
a fixed assignment of subcriticality to the minus channel for
every harmonic. Numeric pressure roots are isolated in the stated
critical neighborhood; high-precision derivatives use mpmath.
See the Python script and its JSON/CSV/plot outputs.

## 3. Extension to exact neutral-plate three-dimensional capillarity

When horizontal plates meet an axially extruded wall at **90° Young
angle**, each regular 2D circular interface extrudes to a true
constant-mean-curvature cylinder satisfying the **three-dimensional**
Young–Laplace equilibrium equation and all contact angles.

At the chosen local pitchfork point, the derivative of the circular
interface family is a positive capillary Jacobi zero mode
\[
\psi(y)=A_M\sqrt{1-y^2/R^2},\qquad
A_M=1+\frac{w(s_*)\,w''(s_*)}
{4[1+(w'(s_*)/2)^2]}.
\]

For \(M\to\infty\), \(A_M\to A_0=3.0957134>0\);
the computed pairs at M=15,51,201 all have \(A_M>3.07\).
The zero mode is nodeless and thus the lowest planar eigenmode
under the free-contact-line Robin conditions. The vertical
Neumann modes add \((n\pi/d)^2\) and **cannot destabilize the
selected critical cylinder at that point**.

Since z-invariant equilibria are an exact invariant subproblem,
and the other depth-varying modes are noncritical, their
local center-manifold elimination preserves the planar
quartic sign at neutral plates, scaled by depth. For the
sampled M values, the analytically constructed circular interfaces
also pass a sampled nonlocal wall-clearance check.

**Boundary of this result:** This is a *special exact stationary
3D geometry with a local spectral argument*, not a 3D
non-neutral-plate numerical solution. An earlier instability
elsewhere on the path, non-graph interfaces, actual hysteretic
chemical pinning, and finite snap landing are still unresolved.
At the manufactured scales for large \(M\), continuum wetting
and smooth-wall assumptions are physically meaningless.

## 4. Independent 3D quartic extractor: implemented but NOT validated

A separate file, `experimental_quartic_kkt.py`, implements the
finite-dimensional constrained fourth-order coefficient using a
bordered KKT system. The Lagrangian stationarity equations are

\[
R(x,p,\mu,c)=(E_x-pV_x-\mu C_x,\ V-V_0,\ C-c)=0.
\]

The envelope theorem gives \(E_{\rm eff}'(c)=\mu(c)\) and hence
\(B_4=\mu'''(0)/24\). Automatic differentiation and
three bordered linear solves compute \(\mu'''\) without fitting
noisy \(c^4\) energy differences. An **independent synthetic
test with an exactly known Schur-reduced quartic coefficient**
passes to double precision.

The exploratory adapter
`experimental_3d_quartic_diagnostic.py` connects this
extractor to the full nonlinear surface energy of Round 6
for two reflected menisci sharing the same total liquid volume.
The work is **not yet a validated prediction**: local tests
show that shallow sidewall contact angles cause underresolved
contact geometry, that stationary roots can shift branches,
and that raw coarse-grid quartic signs can change under
mesh modification even when the nonlinear KKT residual is tiny.

This is an explicit **failed mesh-convergence gate**, not
a result to conceal. No numerical non-neutral 3D \(B_4\)
values are endorsed. Required next tests:

- Track the local symmetric branch in *volume* and solve
  \(B_2(V)=0\) for **each mesh**.
- Resolve near-wall contact-line cusping with clustered /
  adaptive elements and higher-order geometry representation.
- Project out gauge / fixed-volume directions; certify all
  orthogonal constrained eigenmodes and KKT conditioning.
- Verify the **90° exact continuum 3D quartic**
  against \(d B_4^{2D}\), then repeat at 89° and 88°,
  including multiple depths and mesh orders.
- Cross-check in an independent FEM/free-contact-line framework.

The heavy 3D JAX prototype is not called by the routine
GitHub Actions workflow. It is a reproducible experiment,
**not a passed benchmark**.

## 5. Relation to existing literature and evaluation

- [Grimm & Baake, Homometric Point Sets and Inverse Problems (2008)](https://arxiv.org/abs/0808.0094):
  long-established spectral non-identifiability.
- [Roughness contact-line physics, Iliev et al., PRE (2018)](https://doi.org/10.1103/PhysRevE.98.060801):
  nontrivial dependence of contact-line roughness on substrate
  and contact angle, not our particular deterministic example.
- [Wells et al. (2018)](https://doi.org/10.1038/s41467-018-03840-6):
  smooth-topography snapping and pitchfork physics.
- [Bisquert et al. (2024)](https://doi.org/10.1021/acs.langmuir.4c02908):
  competing capillary bifurcations.
- [Original Hatipogullari et al., PRF (2019)](https://doi.org/10.1103/PhysRevFluids.4.044008):
  original capillary statics.

### Professor verdict, round 8

Mathematical theorem in a *specified idealized model*: **8/10**.
Novelty vs inverse-problem mathematics: **moderate**, because
general homometry/phase retrieval is old. New coupling to a
capillary bifurcation is potentially distinctive.
Nonneutral 3D quartic numerics: **NOT VALIDATED**.
Experimental correlation: **NOT VALIDATED**.
Journal submission readiness: **4/10** for capillary physics;
possibly **5/10** for a standalone carefully situated
mathematical note after independent proof review, priority audit,
and verification of all regularity assumptions.

A publication should **never** claim that nanometre/atomic
roughness engineered at arbitrary correlation order is a
physically testable CAH device. Realistic 0.4–0.5 µm
second/fourth-harmonic designs remain the better test case.

Reproduce the theorem and neutral-plate checks:

```bash
pip install -r cah-statics/requirements.txt
python cah-statics/arbitrary_order_capillary_theorem.py
```

The optional JAX KKT extractor requires
`cah-statics/requirements-3d.txt`; run its self-test using

```bash
python cah-statics/experimental_quartic_kkt.py
```

The full 3D diagnostic remains **experimental, manual-only**.
