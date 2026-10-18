# Round 15 — an exact single-front pressure–volume criterion for two-front onset

**Provenance:** this same transfer identity, including its sixth-order degeneracy, was **already derived and tested** in [the unmerged experimental-volume research branch's theory note](https://github.com/mhatipog/scientific-artifact-adoption/blob/research/cah-phase-volume-transfer-20261018/cah-statics/theory/PRESSURE_VOLUME_TRANSFER.md) and [source](https://github.com/mhatipog/scientific-artifact-adoption/blob/research/cah-phase-volume-transfer-20261018/cah-statics/theory/single_front_transfer.py). **It is not newly discovered here.** The new work in the present theory branch is an independent analytic-jet evaluation on the more resolvable Round-14 geometry and an explicit conservative pressure-measurement sign certificate.

**Status:** exact chain-rule theorem **within the previously specified reduced planar fixed-volume model**, checked independently in both signs of the Round-14 Fourier-phase example. **Not** a theorem of arbitrary 3D free-contact-line stability. No direct experimental pressure–volume curve has been calibrated for Yan (2025); the result does not validate one.

## Proposition (local quartic transfer)

Let \(H\in C^3\), \(H'(s_c)>0\), and \(p\in C^3\) be a one-front equilibrium pressure law. Introduce the strictly monotone single-front volume coordinate \(Q=H(s)\), and denote the same pressure as \(P(Q)=p(H^{-1}(Q))\). Assume a symmetric two-front plug obeys

\[
V/d=H(s+c)+H(s-c),
\]
has equal pressures at both fronts, and has the **already established** fixed-volume quartic coefficient

\[
B_4=\frac14\left[H''(s)p''(s)-\frac{H'(s)p'''(s)}{3}\right]_{s_c}
\]
when \(p'(s_c)=0\). The coefficient implicitly uses the capillary-pressure/energy sign convention of the original reduction.

**Then**
\[
\boxed{B_4=-\frac{[H'(s_c)]^4}{12}P'''(Q_c)},\qquad Q_c=H(s_c).
\]
In particular, if \(P'''(Q_c)\ne0\), the sign of the *single-front third pressure–volume derivative* determines the reduced local quartic class: \(P'''(Q_c)<0\) gives \(B_4>0\) (locally supercritical, other modes stable), and \(P'''(Q_c)>0\) gives \(B_4<0\) (locally subcritical, finite-jump-capable but no guaranteed global snap).

**Proof.** Put \(A=H'(s_c)>0\), \(B=H''(s_c)\). Since \(p=P\circ H\) and \(p'=0\), also \(P'=0\), and direct differentiation yields
\[
p''=P''A^2,\qquad p'''=P'''A^3+3P''AB.
\]
Substituting into the specified quartic expression,
\[
B_4=\frac{1}{4}\left[BA^2P''-\frac{A}{3}(A^3P'''+3ABP'')\right]
=-\frac{A^4}{12}P'''.
\]
The \(H''\) contributions cancel **exactly**; no approximation, special Fourier profile or numerical fit is used. \(\square\)

If additionally \(P\in C^5\) near the critical state, an equivalent pressure-balance normal form is obtained by setting \(q=(Q_1+Q_2)/2\), \(\eta=(Q_1-Q_2)/2\). The exact equilibrium condition is \(P(q+\eta)=P(q-\eta)\). Near \(q=Q_c\):
\[
0=2P'(q)\eta+\frac{P'''(q)}3\eta^3+O(\eta^5).
\]
If additionally \(P''(Q_c)P'''(Q_c)\ne0\), the nonzero branch satisfies
\[
\eta^2=-6\frac{P''(Q_c)}{P'''(Q_c)}(q-Q_c)
+O((q-Q_c)^2).
\]
This predicts the *local* branch side, not the finite snap distance or 3D global stability. When \(P'''(Q_c)=0\), the quartic vanishes and higher-order terms matter; the simple sign classification is **inapplicable**, not automatically a transition prediction.

## Numerical audit: a geometrically homometric local pair

For the Round-14 5th/7th phase geometry at *assumed, unmeasured* Young sidewall angle 32°,
\[
w_\pm(x)=1+0.05\cos(2\pi x/0.3)
\pm0.00017[\cos(10\pi x/0.3)-\cos(14\pi x/0.3)]\quad\mathrm{mm},
\]
the exact analytic-jets implementation provides:

| Case | \(s_c\) (mm) | \(H'(s_c)\) | \(P''(Q_c)\) | \(P'''(Q_c)\) | \(B_4\) from both formulas |
|---|---:|---:|---:|---:|---:|
| − phase | 0.7039537629 | 4.035332058 | −13.93848776 | +109.2217413 | −2413.486981 |
| + phase | 0.7088690566 | 2.999132549 | −12.75881297 | −534.8302698 | +3605.930678 |

The older [PR-14 energy-based implementation](https://github.com/mhatipog/scientific-artifact-adoption/blob/research/cah-phase-volume-transfer-20261018/cah-statics/theory/single_front_transfer.py) was additionally rerun: its fixed-volume quartic and sixth-order-degeneracy checks passed, and its coefficient matches the new implementation **exactly to printed double precision for both phase geometries**. Both independent expressions for \(B_4\) agree to better than \(10^{-9}\) at those isolated local roots (observed absolute differences ~\(10^{-12}\)). This is an algebraic cross-check of a specific numerical reduced-model calculation, **not** an independent experimental or continuum three-dimensional prediction.

## A conditional four-pressure-point sign certificate

The theorem suggests measuring \(P'''(Q_c)\) from an independently **calibrated** single-front \(p(Q)\) curve. Because numerical third differentiation amplifies noise, do not classify signs from a smoothing spline without a declared resolution test.

At equally spaced calibrated volumes \(Q_c\pm h,Q_c\pm2h\), define
\[
D_{3,h}=\frac{P(Q_c+2h)-2P(Q_c+h)+2P(Q_c-h)-P(Q_c-2h)}{2h^3}.
\]
If each pressure reading has a *deterministic absolute error bound* \(\epsilon\), and \(P\in C^5\) with an independent **certified** bound \(|P^{(5)}|\le M_5\) throughout \([Q_c-2h,Q_c+2h]\), then Taylor's formula through order 4 and the triangle inequality yield the safe enclosure
\[
\boxed{|D_{3,h}-P'''(Q_c)|\le 3\epsilon/h^3+(17/60)M_5h^2.}
\]
The \(17/60\) is deliberately conservative: the four fifth-order remainders contribute \((32+2+2+32)M_5h^5/(120\cdot2h^3)\). Therefore only \(|D_{3,h}|\) **strictly above** this error budget certifies the derivative sign, and the opposite quartic sign. An interval that overlaps zero must return **indeterminate**. Volume-coordinate error, correlation of errors, uncertain location of \(Q_c\), incomplete equilibration and wetting hysteresis must also be bounded independently; the displayed criterion does **not** address them.

The code demonstrates this bound on synthetic cubic–quintic pressure curves only. **No current Yan pressure/image pair passes the prerequisites**: volumetric liquid thickness is not calibrated, the camera/pressure clock match is unverified, 95° nominal wettability is not independently measured advancing/receding Young angle, and a finite \(M_5\) is not available.

## Context, novelty and gates

Wells et al., *Nature Communications* **9**, 1380 (2018), doi:10.1038/s41467-018-03840-6, and Bisquert et al., *Langmuir* **40**, 24387 (2024), doi:10.1021/acs.langmuir.4c02908, already establish the importance of snaps and competing pitchfork/fold transitions. Therefore neither symmetry breaking nor smooth-geometry hysteresis is claimed here as a discovery. The chain-rule cancellation is a concise **model-specific reparameterization theorem already present in our own draft PR #14**, not a fresh result of this round; a comprehensive mathematical priority search is still required before presenting it as novel to the literature. Exact higher-order geometric correlations of Round 14 are a separate hypothesis of practical interest.

1. Reproduce the unit tests and both numerical signs, offline.
2. Validate the local fixed-volume Hessian and quartic against a converged 3D Young–Laplace free-contact-line calculation for neutral plates, then non-neutral plates. Explicitly test depth modes; a finite polynomial basis is not continuum certification.
3. Measure a single-front \(P(Q)\) with independent volume, pressure, geometry, advancing/receding angles and synchronization; obtain a credible uncertainty bound on its *third derivative*.
4. Use that input **blindly** to predict the two-front class in a fabricated matched-phase pair; only then assess independent experimental novelty.

**Research verdict:** mathematical rigor 8/10 for the reduced identity and algebra; relative novelty 5/10 pending a serious prior-art search; experimental validation 4/10 for the programme but 0/10 for this bridge; scientific impact potential 7/10; current journal-submission readiness **5/10, unchanged**.
