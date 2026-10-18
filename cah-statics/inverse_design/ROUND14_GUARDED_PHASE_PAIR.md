# Round 14 — Guarded isocorrelated geometry pair and competing-fold limit

**18 October 2026. Scope:** exact low-order Fourier correlation identity; independently implemented and tested *local planar* Young–Laplace normal form; exact stationary neutral-plate 3D cylinders with sampled local Jacobi and wall-clearance checks. **No fabricated geometry, measured contact angle, non-neutral 3D convergence, global snap prediction, or experimental validation of onset class.** Submission readiness remains **5/10**.

## A different question from Round 13

[Round 13](https://github.com/mhatipog/scientific-artifact-adoption/blob/research/cah-phase-volume-transfer-20261018/cah-statics/experimental_volume/ROUND13_CRITICAL_DEPTH.md) shows that the original Yan optical projections cannot yield a physical volume without calibrating *incremental* water thickness. That is a single-front experimental gate; it must never be fed directly into a two-front pitchfork model.

We therefore revisited a separately testable piece of the existing theory: can [Round 7](../ROUND7_HIGH_ORDER_ISOSPECTRAL.md)'s exact fifth-order-isocorrelated *two-front* geometry be enlarged for experimental discrimination **without introducing competing capillary folds?** Simply increasing its high-harmonic wall amplitudes turns out to be unsafe.

## The selected conservative pair

For a reflection-symmetric channel, use millimetres and define

\[
w_\pm(x)=1+0.05\cos(2\pi x/0.3)
 \ \pm 0.00017[\cos(10\pi x/0.3)-\cos(14\pi x/0.3)].
\]

This is **the same 5th/7th harmonic family already proved in Round 7**, but with a carefully selected nominal amplitude **170 nm** rather than 75 nm. It is not a new Fourier theorem. The high modes have axial wavelengths **60** and **42.86 µm**; the maximum **difference in total channel width** between paired profiles is **0.657466 µm**, versus 0.290 µm for the original 75 nm pair. The lateral displacement of a *single wall* relative to the common sinusoidal baseline is no more than **0.1644 µm**.

The two *nominal* profiles have exactly equal:
- mean width and extrema 0.95–1.05 mm, since \(|b|<a/24\) guarantees the fundamental extrema remain global;
- Fourier power spectrum, RMS width and RMS slope;
- **every translation-invariant n-point height-correlation function for n=1,2,3,4,5**, for arbitrary offset vectors.

The exact proof is discrete: among signed Fourier modes \(\pm1,\pm5,\pm7\), no zero-sum tuple of length at most five has an odd number of 5/7 modes. Thus flipping the sign of both high-mode coefficients does not change any Fourier product contributing to those correlations. Such a tuple **does** first appear at sixth order, e.g. \((1,1,1,1,1,-5)\). The exact combinatorial parity condition is independently tested by integer-mode enumeration, without sampled-grid estimates. Sign-reversed profiles are not translations/reflections, since preserving the positive fundamental cosine fixes the translation phase.

**Crucial distinction:** Fabrication errors that independently change the 5/7 amplitudes in the two channels destroy *exact* equality of their high-mode power and correlations. Our separate tolerance scan tests persistence of the **capillary signs**, not persistence of perfect real-world statistical matching. Profile metrology, not a nominal CAD file, must establish actual geometric similarity.

## Local capillary calculation and numerical values

We independently reimplemented the previously specified planar circular-interface pressure and fixed-volume primitive:

\[
p(x)=\frac{2\cos\{\theta_Y+\arctan[w'(x)/2]\}}{w(x)},
\quad
H(x)=\int^x w(\xi)\,d\xi+
w(x)^2\frac{2\alpha-\pi+\sin(2\alpha)}{8\cos^2\alpha},
\quad \alpha=\theta_Y+\arctan(w'/2).
\]

At the isolated local pressure-stationary root \(p'(s_*)=0\), evaluate

\[
B_4=\frac14\left[H''p''-\tfrac13H'p'''\right]_{s_*}.
\]

The implementation evaluates \(w,w',\ldots,w''''\), \(p',p'',p'''\) and \(H',H''\) **analytically**, not by noisy fourth-order differences of a fitted energy. As independent checks, finite-differencing the analytic lower derivatives yields relative errors below \(2\times10^{-6}\), and the earlier 75 nm / 36° Round-7 case is reproduced: \(B_{4,-}\approx-1505.931\) and \(B_{4,+}\approx812.389\).

For the **assumed, not experimentally measured** intrinsic sidewall Young angle \(\theta_Y=32^\circ\):

| Local quantity | Minus phase | Plus phase |
|---|---:|---:|
| Unique stationary root in predeclared critical window | 703.953763 µm | 708.869057 µm |
| \(p''(s_*)\) | −226.9730 | −114.7629 |
| \(B_4\) | **−2413.487** | **+3605.931** |
| Local normal-form class, given other modes stable | **Subcritical** | **Supercritical** |
| Neutral-plate planar Jacobi mode factor | 4.1543 | 3.0994 |

The predicted difference between critical half-lengths is **4.9153 µm**. This is an exact output of the specified reduced model, *not* an observed displacement or pressure extremum. It could motivate a future blinded position readout at submicrometre precision. Subcritical \(B_4\) makes finite jumps **possible**, not inevitable; competing folds, path stability, hysteretic contact-line pinning, and finite landing states remain unresolved.

## Why "increase the corrugation" fails: fold creation

The local pressure stationary points themselves obey \(p'(s;b)=0\). A pair of extra stationary points is born where

\[
\boxed{p'(s;b)=p''(s;b)=0.}
\]

This is an ordinary saddle-node of the **stationary-root equation** and is distinct from the desired volume-controlled pitchfork. At \(\theta_Y=34^\circ\), \(a=0.049\) mm, period 0.303 mm—the most adverse sampled manufacturing corner—we directly solved the two equations for the positive phase:

\[
s=0.702834297\ \mathrm{mm},\qquad
b_{\rm fold}=0.241991386\ \mathrm{\mu m},
\]

with combined residual \(<1.5\times10^{-11}\) in raw equation units. Immediately below it, one stationary root exists in the declared 90-µm interval; above it, three are present. This **explains an actual failure mode**, rather than assuming harmonic signals can be strengthened without consequence.

### Stress-grid design decision

For the pair, predeclare the following limited uncertainty ranges:

- Sidewall intrinsic Young angle **30°, 31°, 32°, 33°, 34°**, *not* measured;
- fundamental amplitude \(a=0.049,0.050,0.051\) mm (±2%);
- period \(0.297,0.300,0.303\) mm (±1%);
- **independent** absolute high-mode amplitude errors \(0.8b_0,b_0,1.2b_0\) for each of the two phase profiles (±20%).

Exactly **405 joint parameter combinations** were evaluated, each by isolating stationary roots on \([0.65,0.74]\) mm scaled with the period. At the **selected** nominal \(b_0=0.170\) µm:

| Gate, within the discrete specified grid | Result |
|---|---:|
| Grid points preserving one isolated local root in each profile and opposite B4 signs | **405 / 405** |
| Smallest signed quartic margin \(\min(-B_{4,-},B_{4,+})\) | **239.869** |
| Smallest magnitude \(-p''\) | **43.549** |
| Smallest neutral-plate Jacobi factor | **2.742** |
| Sampled circular-interface wall clearance (excluding wall-contact node) | **positive at all 405 combinations** |
| Worst allowed positive-phase amplitude \(1.2b_0\) | **0.204 µm** |
| Margin to computed worst-corner new-root threshold | **37.991 nm** |

The wall clearance number depends on the interior sampling mesh and **must not** be interpreted as a finite continuum clearance margin: clearance approaches zero at the prescribed contact point. The grid itself is **not an interval enclosure** over the full continuous uncertainty box.

An aggressively enlarged nominal \(b_0=0.200\) µm passed the same **405 sampled** sign/root checks, but its worst positive amplitude, 0.240 µm, lies just **1.991 nm below** the fold-creation point. At \(b_0=0.210\) µm, **6 of 405** combinations already fail the single-root gate. The selected 170 nm candidate deliberately sacrifices some absolute wall-profile contrast to retain approximately **38 nm** of modeled fold headroom. This is evidence of a nontrivial **signal-versus-root-proliferation trade-off**.

[Deterministic summary](derived/guarded_pair_summary.json), [45-stratum aggregate across the 405 combinations](derived/guarded_pair_45_strata.csv), [vector profile/bifurcation figure](figures/phase_pair_robust_window.svg) and [offline SVG generator](plot_phase_pair_window.py) are provided for review.

## Exact neutral-plate statement versus missing finite-depth validation

With plate Young angle **90°**, the circular planar menisci extrude into exact stationary 3D cylinders, and normal perturbations separate into vertical Neumann modes. The minimum computed nodeless-Jacobi factor in the selected stress grid is \(2.742>0\). Under the same Sturm–Liouville and other-mode assumptions explicitly described in [Round 8](../ROUND8_ARBITRARY_ORDER_THEOREM.md), depth-dependent modes are noncritical at the chosen cylinder. This supports the **special neutral-plate local lift**, *not* a general finite-depth non-neutral result or global stability of the actual liquid plug.

A non-neutral top/bottom plate bends the free contact line and invalidates the exact cylindrical solution. Prior direct nonlinear 3D quartic computations failed mesh-convergence checks; this round **does not repair them**. The true experimental bifurcation class remains unvalidated. Also missing are a complete fixed-volume four-degree-of-freedom Hessian audit away from the center mode, branch tracking through competing folds, independent microscopic wetting measurements, chemical heterogeneity, fabrication uncertainty metrology, and a blinded experiment.

## Reproduction, code integrity, and prior art

Run **offline** (requires NumPy and SciPy, no original videos or Actions):

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_phase_pair_window.py' -v
python cah-statics/inverse_design/phase_pair_window.py
python cah-statics/inverse_design/phase_pair_window.py --grid-csv > /tmp/cah_pair_grid.csv
python cah-statics/inverse_design/plot_phase_pair_window.py
```

[The complete source](phase_pair_window.py) and [eight network-free tests](tests/test_phase_pair_window.py) verify exact mode parity, the older 75 nm independent benchmark, analytic-vs-finite-difference derivative agreement, extrema, nominal sign/Jacobi/clearance, all 405 discrete parameter combinations, the computed competing fold, and valid reproducible SVG XML. Test output: **8/8 PASS**; vector XML parsed successfully. No workflow changes, workflows manually dispatched, original microscope media, or extra dependencies were added.

**Prior art / originality boundaries:**
- [Grimm & Baake (2008), *Homometric Point Sets and Inverse Problems*](https://arxiv.org/abs/0808.0094): spectral-statistical nonuniqueness is long established. The fifth-order profile construction itself appears in **our own Round 7**, not in this round for the first time.
- [Wells et al. (2018), *Nature Communications*](https://doi.org/10.1038/s41467-018-03840-6) already describe sinusoidal-topography snaps and subcritical pitchforks.
- [Bisquert et al. (2024), *Langmuir*](https://doi.org/10.1021/acs.langmuir.4c02908) explicitly identify competing folds/pitchforks controlling symmetry breaking.
- [Kumar & Harvie (2023), *Contact Angle Hysteresis on Rough Surfaces Part II*](https://arxiv.org/abs/2303.09149): micro-roughness energy dissipation and open-source Surface Evolver modeling give a strong established baseline for future full 3D comparisons.

Our specific contribution **this round** is *an independently coded, conservative geometry-design window* exposing a fold-proliferation trade-off and a candidate that is more optically distinguishable, **not a new theorem of homometry or the discovery of capillary snaps**. The literature search did not identify this exact phase-pair criterion, but absence in a search is not a priority proof.

## Professor's assessment — Round 14

| Category | Score | Specific justification |
|---|---:|---|
| Mathematical rigor | **8/10** | Exact fifth-order Fourier parity from Round 7 and independent derivative and root-fold calculations; no continuous interval certificate or complete constrained Hessian |
| Novelty relative to existing literature | **5/10** | More realistic and falsifiable design trade-off, but standard homometry/fold/pitchfork mathematics and no proven first-in-field result |
| Experimental validation | **4/10** | Existing original Yan-image experiments remain useful but are unrelated to this two-front pair; **zero actual validation** of this new design |
| Potential impact | **7/10** | A metrology-scale paired-channel experiment could distinguish whether high-order Fourier phase predicts capillary onset |
| **Current journal-submission readiness** | **5/10 unchanged** | No fabricated-channel test or converged non-neutral finite-depth pressure/onset result |

**Strongest reviewer objection:** This is a continuum idealization of two extraordinarily similar smooth walls; fabrication/chemical disorder may exceed their nominal difference and a local quartic sign cannot establish the observed global hysteresis. The ±20% amplitude study destroys exact fifth-order spectral matching when errors differ, further requiring measured-profile verification.

**Single most decisive next calculation:** Perform *mesh-, order-, and fixed-volume-branch-converged* free-contact-line 3D Young–Laplace continuation on **the exact selected 170 nm pair**, first reproducing the analytic 90°-plate quartic **and then at 89° and 88°**, with independently specified sidewall Young angles and constrained Hessian eigenvalue checks. If its opposite onset classes persist, preregister a blinded matched-PSD channel fabrication/imaging experiment. **Potential journal:** *Langmuir* or *Physics of Fluids* after that gate; *Physical Review Fluids* requires an independently validated general principle. Do not change the 5/10 submission readiness just because the grid grew.
