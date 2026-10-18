# Round 13 — A depth-threshold certificate from original 95° frames

**18 October 2026. Status:** New *conditional, experimentally grounded conservation bound*; not a measured 3D volume, not an independently validated capillary-pressure prediction, and not a novel observation of Haines jumps. Publication readiness **5/10, unchanged**.

## Question and source provenance

We need a test that can be completed **before** solving the inverse 3D wetting geometry and without assigning an arbitrary camera–pressure time shift. Use the original Yan et al. microscope TIFF (DaRUS data file [351440](https://darus.uni-stuttgart.de/api/access/datafile/351440), [doi:10.18419/DARUS-4683](https://doi.org/10.18419/DARUS-4683), CC BY 4.0; image-derived quantities only). The TIFF archive name says 90° but its embedded Nikon *source folder* says `Theta_95/Q1`; see [Round 12](ROUND12_SPARSE_95_CLOCK.md). The archive MD5 is **not independently verified from a whole-file download**. The 95° camera/pressure same-run clock and acquisition-date discrepancy remain unresolved; **no pressure trace is paired**.

The nominal pump rate is `0.01 mL/h = 0.0027777778 µL/s`. The image pitch is `2.4859 µm/px`, and the nominal maximum liquid column depth is `D = 120 µm`. Original pixel pore walls are `(116,69)–(438,391)`; three-pixel wall strips are excluded, yielding a `316 × 316` analysis ROI. Five thresholds are **frozen from Round 12** and were not tuned to make the new depth threshold small. Optical darkening is compared pixelwise with original first and last frames after Gaussian blur σ=1 and static exterior drift correction.

The independent sparse byte-range audit in [`critical_depth_certificate.py`](critical_depth_certificate.py) requests strictly verified HTTP 206 ranges for the TIFF header, first/final frame and burst frames 4818–4824 (**approximately 4.72 MB**, not the full 2.54 GB). The result is in [`derived/critical_depth_95_windows.csv`](derived/critical_depth_95_windows.csv); no original TIFF is republished and no Actions workflow is triggered.

## A short conditional theorem: activation–loss conservation

Consider two registered projected masks, separated by `Δt`. Let `A₊` be the physical projected area of pixels changing **dry/negative → wet/positive** and `A₋` the area changing **wet/positive → dry/negative**. These are measured **separately**; net darkening can conceal substantial threshold reversals. Let `Q Δt` be nominal delivered liquid, assuming negligible net losses in transit from the pump control volume. Let `D` bound the maximum water-column depth removable on each disappearing pixel.

Assume that an independent depth-sensitive measurement demonstrates a **lower bound** `h₀` on the **area-averaged increment** in actual liquid-column depth across the newly positive pixels. This is stronger than demonstrating a final (rather than changed) layer thickness: a pixel previously below the optical threshold can already contain liquid. Define `R` as the net additional supply/release into the newly positive area from *all remaining* compartments, including thinning of pixels whose class stays positive, subthreshold films, pre-existing throat liquid, and upstream tubing/bubbles; neither its physical location nor compliance law is presumed. Nonnegative outgoing flows only increase the supply demand. Then mass conservation gives

`h₀ A₊ ≤ Q Δt + D A₋ + R`, hence

\[
\boxed{R\ \ge\ \max\{0,\;h_0 A_+ - D A_- - Q\Delta t\}.}
\]

Equivalently, for the **pump plus fully drained disappearing pixels to suffice without supply from any other area**, a necessary condition is

\[
\boxed{h_0\ \le\ h_{\mathrm{crit}}^{\mathrm{guard}}
=\frac{Q\Delta t+D A_-}{A_+}.}
\]

The simpler `QΔt/(A₊−A₋)` assumes that the areal net change is representative of a single depth; it ignores the possibility that disappearing pixels were full-depth and is **not** the rigorous guard used in conclusions. This bound is a conservation inequality, **not an original physical law or a fitted compliance model**. The unmeasured assumption is precisely `h₀`, the incremental water depth per newly positive pixel. If `h₀` cannot be calibrated, `R` is not experimentally lower-bounded above zero. Measuring a final water-film thickness alone is insufficient if the pixel initially contained water.

### Measured classification counts and conditional depth thresholds

Each window contains **two camera steps, Δt = 0.066657007 s**, so nominal pump delivery is only **0.0001851584 µL**. All values below were independently reproduced from original TIFF byte ranges. Optical segmentation count/fraction labels are fixed from Round 12.

| Window | Classifier | Newly positive pixels | Reversed pixels | Reversals/activations | Fully-loss-guarded `hcrit` |
|---|---|---:|---:|---:|---:|
| 4820→4822 | 100 / 0.4 | 27,089 | 134 | 0.49% | 1.700 µm |
| 4820→4822 | 180 / 0.4 | 26,349 | 60 | 0.23% | 1.410 µm |
| 4820→4822 | 180 / 0.5 | 26,945 | 1,153 | 4.28% | 6.247 µm |
| 4820→4822 | 180 / 0.6 | 28,055 | 2,936 | 10.47% | 13.626 µm |
| 4820→4822 | 260 / 0.6 | 27,290 | 2,904 | 10.64% | 13.867 µm |
| **4822→4824** | **100 / 0.4** | **24,580** | **37** | **0.151%** | **1.400 µm** |
| **4822→4824** | **180 / 0.4** | **23,686** | **22** | **0.093%** | **1.376 µm** |
| **4822→4824** | **180 / 0.5 (reference)** | **28,528** | **29** | **0.102%** | **1.172 µm** |
| **4822→4824** | **180 / 0.6** | **33,773** | **38** | **0.113%** | **1.022 µm** |
| **4822→4824** | **260 / 0.6** | **32,930** | **20** | **0.061%** | **0.983 µm** |

The *early* interval is significantly **segmentation-sensitive**: as many as 2,936 pixels reverse classification at strict thresholds. We **reject** treating its small unguarded 1.1–1.2 µm depth as a robust certificate. The *later* interval is more coherent: at most 38 reversed pixels among ≥23,686 newly positive pixels; even allowing **every reversed pixel to drain the entire nominal 120 µm**, its guarded critical increment stays **0.983–1.400 µm** across all five masks.

A deliberately oversized max-depth assumption `D=200 µm` in place of the nominal 120 µm raises the largest later-window guarded threshold only to **1.520 µm**. This probes the loss-pixel credit, *not* previously positive pixels that might thin.

For illustration **only**, if an independent 3D calibration later establishes a lower bound on *area-averaged incremental wet thickness* of `h₀=2 µm` across newly positive pixels, the equation requires `R ≥ 0.000091–0.000207 µL` for the later interval under `D=120 µm`. This is **conditional**, not an inferred or measured stored volume, and in particular is **not evidence that upstream apparatus compliance released this amount**. The data alone permit arbitrarily thin activated films and cannot establish a positive `R`.

[Figure: log-scale threshold comparison](figures/critical_depth_95_guarded.svg), reproducible offline using [`plot_critical_depth.py`](plot_critical_depth.py). The horizontal 2 µm line is marked as a **hypothetical design target**, not an observed pore-film thickness.

## Independent kinetic falsification attempt and why it was withheld

We also explored a completely different, depth- and pressure-clock-free diagnostic: if measured front coordinate `x(t)` itself obeys an **autonomous quadratic** (Riccati) ODE with fixed coefficients, its equally spaced time-step maps are Möbius transformations. Thus four-point position cross-ratios should be time-translation invariant. This is classical projective mathematics, not a new theorem. Using the originally resolved 95° center-band front positions for frames **4812–4824** `[225,230,235,241,247,256,267,281,298,316,333,360,395]` pixels, the step-two cross-ratios include **1.26257 (4816–4822)** and **1.42902 (4818–4824)**. A global least-squares Möbius flow fit nevertheless has **1.45 px RMS residual** and ~2.46 px maximum residual. Early-window extrapolations are highly unstable. Without an independently calibrated <~1 px systematic front-location error, we **cannot claim a robust rejection** even of that directly parameterized normal form.

More fundamentally, front kinematics **cannot identify** a unique capillary pressure–position curve if mobility/forcing are unknown: along any strictly monotone observed `x(t)`, and for almost any proposed `p_eq(x)`, the missing driving pressure or mobility can be chosen to reproduce the trajectory. Fitting a Riccati front is not proof of a saddle-node equilibrium bifurcation. We therefore **withhold that tempting kinetic claim**; this audit prioritizes conservation over overfitted normal forms.

## Comparison with established literature and novelty boundary

- [Yan et al. (2025), *JCIS* 696, 137884](https://doi.org/10.1016/j.jcis.2025.137884) directly report wettability-dependent pore-filling and capillary pressure dynamics, and already speculate that an inlet/sensor air bubble could affect volumetric-flow interpretation. The new result **does not discover** this possibility; it converts a particular observed burst into a *predeclared depth-calibration falsification threshold*, with no pressure time shift.
- [Berg et al. (2013), *PNAS* 110, 3755–3759](https://doi.org/10.1073/pnas.1221373110) already directly imaged fast 3D Haines jumps and energy redistribution. No claim of discovering rapid pore filling, hidden flow or conservation laws is warranted.
- The existing [Round-8 arbitrary-order Fourier-phase theorem](../ROUND8_ARBITRARY_ORDER_THEOREM.md) and the [additive two-front transfer result](../theory/PRESSURE_VOLUME_TRANSFER.md) remain distinct theoretical findings. This result concerns a **single moving interface**, not a two-front pitchfork, and is not a successful 3D Young–Laplace pressure-extremum validation.

This simple certificate improves **falsifiability and experiment design**, not current evidence for a new capillary bifurcation principle. To change that status, depth/pressure must be measured independently, and a pressure prediction predeclared from geometry, independently measured advancing/receding wettability, measured plate separation and rounded throat corner must then survive comparison with the **same-run synchronized** pressure data.

## Reproduction and tests

From the repository root:

```sh
python -m unittest discover -s cah-statics/experimental_volume/tests -p 'test_critical_depth_certificate.py' -v
python cah-statics/experimental_volume/critical_depth_certificate.py --download-ranges > /tmp/critical_depth_check.csv
python cah-statics/experimental_volume/plot_critical_depth.py
```

The first and third commands are offline (the third writes the deterministic SVG); only the **explicit** second command contacts the original public DaRUS archive, with strict per-range HTTP 206 validation and no raw-image persistence. The frozen CSV is checked against an independently executed full range replay. **Seven network-free tests passed**, spanning mask activations and reversals, maximal-depth credit, exterior drift, invalid inputs and a synthetic expanding front. The generated SVG passed XML parsing and the expected ten-bar count. No Actions were dispatched. No original microscope data were committed.

## Professor's assessment

| Category | Score | Evidence |
|---|---:|---|
| Mathematical rigor | **8/10** | Explicit conditional mass-conservation inequality, unit-checked numerical audit, rigorous worst-case disappearing-pixel credit; liquid increment still unmeasured |
| Novelty relative to existing literature | **4/10** | Useful experiment-specific inversion of a standard mass balance, but Haines bursts and bubble reservoirs predate this audit |
| Experimental validation | **4/10** | Original source pixels independently reread, two sign-separated mask windows and five threshold checks; no calibrated z-depth or synchronized pressure |
| Potential impact | **6/10** | Transforms the missing 3D measurement into a practical ~1.5–2 µm falsification gate, portable to other optical pore studies |
| Journal submission readiness | **5/10, unchanged** | Still no novel quantitatively validated pressure–geometry prediction |

**Strongest reviewer objection:** The bound merely reformulates mass conservation and the pivotal `h₀` is not measured. Depending on tracer optics, newly dark pixels might add far less than 1 µm water or might have already contained optically invisible liquid, and front speed says nothing unambiguous about capillary pressure.

**Most decisive next experiment:** On a repeat **95° Q1** run, synchronize pump encoder, high-speed microscope, pressure sensor, and **depth-sensitive local liquid-column measurement**. Precalibrate *incremental* liquid thickness on newly dark columns and thinning on previously dark columns against a static reference, target an uncertainty smaller than **0.5 µm**, and verify an area-averaged increment exceeding the 1.5–2 µm threshold without using pressure to tune segmentation. A calibrated h₀ (and ideally a spatial height map) converts the present conditional certificate into a measured redistribution lower bound. Separately resolve original `test1035.nd2` acquisition time and pairing with the deposited pressure trace before pressure-curve validation. *Langmuir*/*Physics of Fluids* are plausible later targets; *Physical Review Fluids* needs a independently predictive general capillary principle.
