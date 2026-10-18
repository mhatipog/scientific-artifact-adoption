# Round 9 — The singular limit behind arbitrary-order capillary instability

**Date: 18 October 2026. Theory and numerical verification only.**
This is the next step after the arbitrary-order roughness-correlation
construction in Round 8. The new finding changes its *physical*
interpretation: the local quartic sign remains different even as
all reasonable finite-displacement energy signatures disappear.

## 1. Research question

How can two wall profiles converge in \(C^3\), and have identical
all correlations up to any prescribed finite order, while their
local Young–Laplace quartic coefficients remain of opposite sign?

Consider the Round 8 construction at a nondegenerate sinusoidal
quartic degeneracy:
\[
w_{\pm,M}(x)=1+0.05\cos(kx)
 \pm\frac{\beta}{M^4}
  [\cos(Mkx)-\cos((M+2)kx)],
\]
where \(k=2\pi/(0.3\,\mathrm{mm})\), \(\beta=0.046875\) mm,
\(M\) is an odd integer and the intrinsic sidewall Young angle is
the numerically defined sinusoidal degeneracy
\(\theta_*\approx37.145146436^\circ\).

This is a **continuum-model mathematical construction**. As M
increases, the high-harmonic amplitudes become smaller than atomic
length scales. Such geometries are *not* fabricated microchannels.

## 2. Exact fixed-volume energy identity

In the symmetric two-front planar Young–Laplace model define
\[
p(x)=\frac{2\cos[\theta+\arctan(w'(x)/2)]}{w(x)},
\quad H(x)=\int_0^x w(t)dt+
 w(x)^2\,f[\theta+\arctan(w'(x)/2)],
\]
\[
f(t)=\frac{2t-\pi+\sin(2t)}{8\cos^2t}.
\]
Let half-length \(s(c)\) vary while plug-centre \(c\) is displaced
so that exact area (liquid volume divided by channel depth) is fixed:
\[
H[s(c)+c]+H[s(c)-c]=2H(s_*).
\]
The circular-arc capillary-energy and area derivatives satisfy the
exact local thermodynamic identity \(W'(x)=-p(x)H'(x)\).
Consequently, without computing a potentially expensive wall
arclength/elliptic integral, one can integrate the **exact**
constrained energy derivative:
\[
\boxed{
\frac{dE}{dc}=
\{p[s(c)-c]-p[s(c)+c]\}
\frac{2H'[s(c)+c]H'[s(c)-c]}
     {H'[s(c)+c]+H'[s(c)-c]}.
}
\]
We have verified the identity by differentiating the independent
surface-plus-wall energy formula.

At the local pitchfork threshold \(p'(s_*)=0\), expand
\(s(c)=s_*+s_2c^2+s_4c^4+\cdots\).
The exact volume constraint yields
\[
s_2=-\frac{H''}{2H'},
\quad
s_4=-\frac{H''s_2^2+H'''s_2+H''''/12}{2H'}.
\]
Define
\[
a_3=2p''s_2+p'''/3,
\quad
a_5=2p''s_4+p'''s_2^2+p''''s_2/3+p'''''/60.
\]
Then the reduced constrained energy coefficients
\(E(c)=E_0+B_4c^4+B_6c^6+\cdots\) are
\[
\boxed{B_4=-H'a_3/4},\qquad
\boxed{B_6=-\{H'a_5+(3H''s_2+H'''/2)a_3\}/6}.
\]
These give an *independent* derivation of the earlier quartic
criterion and a new compact sixth-order criterion.
All derivatives are evaluated at the selected critical half-length.

The code independently checks the \(B_6\) expansion against
**numerical integration of the exact fixed-volume energy** at tiny
displacements. It does not rely on unstable polynomial fitting.

## 3. The crucial sixth-order divergence

For the high-harmonic construction \(b=\pm\beta/M^4\),
\(w^{(4)}\) has an O(1) change, so \(B_4\) retains an O(1)
difference between the ± channels on nondegenerate harmonic
subsequences. However, \(w^{(6)}\) is O(M²), generating
\[
\boxed{B_6/B_4\sim-(kM)^2/30}
\]
**along harmonic subsequences where the oscillatory leading quartic
factor stays bounded away from zero**.

Define a characteristic crossover length, not a measured slip:
\[
\boxed{\ell_M=\sqrt{|B_4/B_6|}
\sim\sqrt{30}/(kM)
=0.872\,\lambda/M.}
\]

| M | Positive-phase B4 | Positive-phase B6 | Crossover scale |
|---:|---:|---:|---:|
| 15 | +866.502 | −2.76323×10⁶ | 17.71 µm |
| 51 | +506.152 | −2.34371×10⁷ | 4.65 µm |
| 201 | +416.552 | −2.36101×10⁸ | 1.33 µm |

**Important nonuniformity:** When the oscillatory quartic
prefactor is unusually small, subleading contributions can dominate;
for example, M=75 does **not** obey the simple sixth/quartic sign
expectation. The asymptotic is not claimed uniformly over *all*
harmonics.

## 4. Beyond Taylor coefficients: a universal finite-displacement law

Let \(\xi=Mkc\). At fixed \(\xi\), along the same
nondegenerate subsequences, a matched high-frequency expansion of
the pressure-volume identity above gives
\[
\boxed{
E_M\!\left(\frac{\xi}{Mk}\right)-E_M(0)
=\frac{24 B_{4,M}}{(Mk)^4}
\left[\cos\xi-1+\frac{\xi^2}{2}+o(1)\right].
}
\]

The universal crossover kernel
\[
\boxed{\Phi(\xi)=\cos\xi-1+\xi^2/2}
\]
has \(\Phi(\xi)=\xi^4/24-\xi^6/720+\cdots\).
It recovers the quartic and sixth-order coefficients automatically.

The mechanism is a singular balance: the periodic high-harmonic
pressure force contributes a \(\sin\xi\) term, but centering at the
*shifted* pressure extremum subtracts its linear part. The
constrained force is therefore proportional to
\(\xi-\sin\xi\), whose integral is precisely \(\Phi\).

**This is not just an asymptotic guess.** Independently solving the
two-front equal-volume equation at every integration point, the
high-precision script evaluates the exact constrained energy.

| M | Exact / universal energy at ξ=2 | E/(γd), mm at ξ=2 |
|---:|---:|---:|
| 15 | 0.993570 | 1.23846×10⁻⁶ |
| 51 | 0.971126 | 5.29122×10⁻⁹ |
| 201 | 1.005075 | 1.86794×10⁻¹¹ |

At M=201 and ξ=2, the full energy agrees within **0.51%** of the
universal expression. The energy difference is in mm because
E/(surface tension × channel depth) is a length.

The difference in exact finite-displacement energy shrinks as
\(M^{-4}\), even while the *infinitesimal* quartic coefficient
difference remains O(1).

## 5. The correct functional-analytic interpretation

The width perturbations vanish in every \(C^r\) norm for
\(r\le3\): \(\|w_{\pm,M}-w_{\rm sinusoid}\|_{C^r}
=O(M^{r-4})\). Yet \(B_{4,\pm}\) stay different along
appropriate subsequences.

Hence, for smooth wall profiles close to a degenerate
bifurcation, the map "geometry → local quartic coefficient"
**cannot be continuous in the \(C^3\) topology**: its formula
depends on the wall's fourth derivative.

But this does **not** make the physical energy itself unstable
under small displacements. The exact energy contrast uniformly
decays in the rapidly oscillatory construction; it is the
**infinitesimal classification**, on a shrinking displacement
window, that is a singular observable.

This distinction is essential for experimental interpretation:
matching finite-order roughness correlations is insufficient
to determine an ideal infinitesimal normal-form sign, but it
does **not follow** that experimentally resolvable jumps remain
different when the geometric modification is arbitrarily small.

## 6. Corrected 3D quartic numerical adapter

An independent audit found a source-model mismatch in the
experimental 3D bordered-KKT quartic prototype: it previously
seeded critical half-lengths for **5th/7th-harmonic profiles at
36°**, but actually invoked the **2nd/4th-harmonic nonlinear
surface solver at 36.5°**.

That was incorrect and could produce spurious numerical
instability. The prototype now computes the actual model's
planar critical seeds from its own width/Young-angle functions
and checks them against:
- negative phase: **0.700331099512868 mm**
- positive phase: **0.697477189842539 mm**

**This fixes a definite source mismatch, not the deeper 3D
convergence problem.** The discretized critical *volume* may
differ from the continuum root. The non-neutral full 3D B4
still has NOT passed branch selection, mesh convergence and
competing-eigenmode certification. Its raw numbers must not
be used in a manuscript.

## 7. Professor assessment and publication positioning

Prior art already includes [Wells et al. (2018)](https://doi.org/10.1038/s41467-018-03840-6)
and [Bisquert et al. (2024)](https://doi.org/10.1021/acs.langmuir.4c02908)
on topography-driven pitchfork and fold instabilities. General
spectral non-identifiability/homometry also predates our research,
e.g. [Grimm & Baake (2008)](https://arxiv.org/abs/0808.0094).

Potential genuine mathematical contribution:
**a quantified singular limit linking high-order roughness
statistics, capillary normal-form sensitivity and vanishing
finite-amplitude energy observables**. The full energy kernel
is a stronger and more physically honest finding than just
the arbitrary-N quartic sign reversal.

Mathematical value: **8.5/10** for a self-contained model theorem.
Journal novelty certainty: **5/10** (priority not exhaustively
established). Physics/experimental submission readiness:
**4/10**, essentially unchanged. Plausible future targets:
*Physics of Fluids*, *Langmuir*, or *Physical Review Fluids*
after independently converged 3D finite-depth stability and
measured geometry/contact-line trajectories. A well-reviewed
mathematical short note may be possible earlier.

Best immediate next experiment is still the physically
resolvable 2nd/4th-harmonic geometry (±0.4–0.5 µm), not the
subatomic arbitrary-N limiting sequence. Compare finite
measurable energy barriers and actual stick-slip trajectories,
not only formal quartic coefficients.

### Reproduce

```bash
pip install -r cah-statics/requirements.txt
python cah-statics/singular_limit_energy.py
```

Creates `round9_singular_limit.csv`,
`round9_energy_collapse.csv`, JSON results, and two
self-contained model-only figures under `cah-statics/outputs/`.
The heavy JAX 3D prototype remains manual-only.
