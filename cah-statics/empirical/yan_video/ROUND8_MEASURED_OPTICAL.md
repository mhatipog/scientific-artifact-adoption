# Round 8 (18 October 2026): quantitative filling from original microscopy

**Research status: genuine independent raw-image measurement, not a parameter-free
validation of our capillary-pressure or isospectral phase-switch theory.**

## Sources, data integrity, calibration

Original CC-BY 4.0 experimental data:
Müller, Yan, Raoof & Weigand (2025),
https://doi.org/10.18419/DARUS-4683;
associated JCIS paper: https://doi.org/10.1016/j.jcis.2025.137884.

We retrieved two Q1 water-air ImageJ TIFF stacks through GitHub Actions,
validated their published original-file MD5 values, and inspected
ALL images in each stack for optical filling (only 12 sampled frames
are retained as a public contact-sheet artifact per experiment).

| Parameter | 60° experiment | 120° experiment |
|---|---|---|
| Source | 60_Q1_water_air_3330-3710_TD.tif | 120_Q1_water_air_4700-5357_TD.tif |
| Public datafile ID | 351428 | 351438 |
| Original bytes | 199,824,766 | 345,104,065 |
| MD5 verified | 3dbfc99446d61b89c905057cd360d5b7 | 1beeb40ae47b8c7299aba034e40c5c5f |
| Frames (512 × 512 pixels) | 381 | 658 |
| TIFF ImageJ frame interval | 0.03332186 s | 0.03333214 s |
| Published spatial sampling | 2.4859 µm/pixel | 2.4859 µm/pixel |

For each experiment, restrict to the **same interior-pore ROI**
x = 110–404, y = 95–389 in the original 512-pixel image;
the result is **not** yet a calibrated independent wall segmentation.
Let I0 and If be initial and final TIFF intensities. Define an optical
signal from pixels that darken by at least a chosen threshold,

\[
S_{\rm opt}(t)=
\frac{\sum_{i\in\Omega}[I_i(0)-I_i(t)]}
     {\sum_{i\in\Omega}[I_i(0)-I_i(t_{\rm final})]}.
\]

Use three independent threshold ratios (0.12, 0.22, 0.35) relative to
the median positive initial-final change, with a minimum contrast of
40 instrument intensity counts. Smooth over 1.4 frames and take the
first crossing of the monotone envelope at Sopt = 0.10 and 0.90.

The three-threshold range is a **method-sensitivity diagnostic**,
not a complete imaging uncertainty estimate or a confidence interval.

## Results: measured timescale versus nominal pump delivery

| Directly observed or derived quantity | 60° | 120° |
|---|---:|---:|
| 10% optical filling time, clip-relative | 6.489 s | 6.594 s |
| 90% optical filling time, clip-relative | 9.976 s | 7.708 s |
| **10–90% duration** | **3.488 s** | **1.114 s** |
| Three-mask duration range | 3.470–3.518 s | 1.110–1.119 s |
| Reported Q1 pump flow | 0.01 mL/h | 0.01 mL/h |
| Pump delivery during 10–90% | 0.009688 µL | 0.003094 µL |
| Nominal pore volume 800 × 800 × 120 µm³ | 0.0768 µL | 0.0768 µL |
| **0.8 × nominal volume / pump delivery** | **6.34×** | **19.86×** |
| Mask-sensitivity range, same ratio | 6.29–6.37× | 19.76–19.93× |
| Reported final pore saturation from published paper | ≈100% | ≈92% |

The 120° optical filling process is ~3.13 times faster than the 60°
case in this particular 10–90% comparison. The uncorrected nominal
flow-ratio 19.86× at 120° is an *overestimate* if residual corner air
leaves only 92% final saturation. A simple 0.92 multiplier gives
~18.27×, still far greater than the imposed direct pump rate.

Under a **quasi-2D uniform-depth volume approximation** the amount
that has entered the pore from previously stored fluid, not directly
delivered by the pump over this event interval, is approximately

- 60°: 0.8×0.0768 − 0.009688 = **0.05175 µL**.
- 120°: 0.8×0.92×0.0768 − 0.003094 = **0.05343 µL**.

These are **conditional estimates**, *not* directly measured gas
volumes or apparatus compliance. The optical signal is not an
independently verified liquid-volume calibration. Volumetric flow
could also reflect tubing elasticity, redistribution, pumping
transients, or partially wetting corners. The original authors
already discussed trapped gas and pressure/flow mismatches; a fast
filling event is not a novel physical discovery of our model.

The two optical curves and original frame contact sheets are in
archived GitHub Actions artifacts (7-day retention):
- 60°: https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39106042033
- 120°: https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39106329140

## Critical pressure–video synchronization falsifier

The pressure sensor timestamps count from its own start. The TIFFs
preserve original *video frame numbers*, but contain **no common
absolute acquisition clock**. Assuming the same time origin for
both instruments would map the optical 10–90% intervals to:

| Case | Hypothetical absolute camera time | Pressure change on that naive alignment |
|---|---|---|
| 60° | 117.45–120.94 s | ~−0.0085 kPa (nearly flat) |
| 120° | 163.26–164.37 s | ~+0.0171 kPa (small rise) |

In contrast, the 120° pressure recording has a prominent
~−0.36 kPa five-second decline at sensor ~165–170 s, and the
largest five-second decline in the 60° recording is near
sensor ~879–884 s. The superficial overlap for 120° and
lack thereof for 60° **do not justify inventing a common offset**.

**Result:** The two streams cannot yet be matched event-by-event
without independent camera–pressure time registration.
The published manuscript describes synchronization but the
public TIFF files examined here lack an absolute camera clock.
Do not label pressure drops as measured 10–90% filling pressure
until clocks are reconciled using acquisition-system records or
predeclared identifiable landmarks.

## What this validates and what it does not

**Validated:** both original optical stacks can be processed
reproducibly, and pore filling is much faster than the nominal
pump-only volume budget under the quasi-2D interpretation; the
timing contrast across wettabilities survives three optical
thresholds.

**NOT validated:** any 3D capillary pressure or critical pressure
prediction, our two-front symmetry-breaking normal form, or the
isospectral phase-dependent classification. Yan et al.'s setup is
a **single-front throat→pore** expansion. It is a worthwhile control
for a future adapted model, not a direct test of the two-front
pitchfork. Event speeds also require hydraulic resistance and
apparatus storage, beyond quasistatic capillary energy.

**Professor assessment:** internal reduced-model mathematics ~7/10;
raw experimental reproducibility ~6/10; parameter-free burst-pressure
validation ~2/10; stand-alone journal-submission readiness **5/10,
unchanged**. No novelty claim for volume-storage physics: authors
recognized it in 2025.

**Next strict gate:** independently infer pore inlet geometry/contact
angles; register camera and sensor clocks; estimate actual pore liquid
volume from segmented phase masks and 3D height; then predict
critical pressure *prior* to looking at aligned measured events and
compare with the authors' existing level-set / Crisp–Thorpe baseline.
