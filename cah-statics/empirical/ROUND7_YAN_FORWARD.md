# Round 7 — Yan single-front adaptation: measured inputs versus assumptions

Status: **partial quantitative forward model / empirical alignment incomplete**.

Original experiment: Yan, Müller et al. (2025), J. Colloid Interface Sci.
[doi:10.1016/j.jcis.2025.137884](https://doi.org/10.1016/j.jcis.2025.137884).
Data archive:
[DaRUS doi:10.18419/DARUS-4683](https://doi.org/10.18419/DARUS-4683).

## Independent source geometry and experimental pressure

Paper Fig. 1 and Sec. 2.1: a 150 µm inlet/outlet connected to
an 800×800 µm square pore, 120 µm depth, entrance corner radius
roughly 20 µm by SEM, wall irregularities ±5 µm. Imaging ~2.48 µm/px
and frames at 30 fps. These are source-derived values, not fitted.

A separate passage in Sec. 4.2 refers to a characteristic rounded
"pinning" region radius **Rc≈5.6 µm**, distinct from the SEM caption
**≈20 µm**. Both occur in the published text. This ambiguity
prevents claiming a uniquely known local wall profile at the contact
line; 5.6 may reflect a separate effective definition or measurement.
Do NOT quietly replace the actual corner with one chosen radius.

Original pressure files, saved from the public DaRUS experiment and
read by the instrument-aware CI parser, contain >125,000 valid samples
across the 60,95,120° low-flow water-air trials. Unit is **kPa**;
AlgorithmUnit=psi refers to a sensor calibration algorithm.

Measured one-second-averaged *candidate* abrupt decreases:
- 95°: 0.94639→0.58476 kPa = **0.36163 kPa**,
  from t≈256.42 to 258.41 seconds after the pressure logging start.
- 120°: 1.01423→0.62797 kPa = **0.38626 kPa**,
  from t≈166.15 to 169.13 seconds after start.
- 60°: multiple unrelated-looking large transients; **no unique
  image-correlated event selected**.

These are pressure features, **not** confirmed pore-filling bursts.
Reproducibility and event uncertainty must be established after image
alignment. The pressure sensor has ±200 Pa total reproducibility
uncertainty; do not confuse a 13 Pa sensor resolution with accuracy.

## One-front nominal rounded-wall baseline

Use an idealized 20-µm quarter-circle fillet at each side of the
150-µm channel:

\[
w(\phi)=150+2R(1-\cos\phi),\quad
x=R\sin\phi,\quad 0\le\phi\le90^\circ.
\]

Widths and lengths above in micrometres. Construct a **signed** static
capillary-pressure profile with \(\gamma\approx72\) mN/m:

\[
P_w-P_a=-2\gamma\left[
\frac{\cos(\theta_\mathrm{side}+\phi)}{w(\phi)}
+\frac{\cos\theta_\mathrm{plate}}{120\ \mu\mathrm m}
\right].
\]

The term in channel depth is a **rectangular-cross-section curvature
closure**. This is **not** the full 3D constant-mean-curvature solver
and does not cover a meniscus entering the 800 µm pore.

The angle can be varied separately on side and horizontal walls.
**Important:** The paper measures *static sessile-drop angles*
60/95/120°, but uses *effective static advancing angles* for low
flow in its simulations. These are not interchangeable. In our script
the side and plate angles equal the stated static angle *only as an
explicit hypothesis*. No pressure data are used to infer them.

| Static-angle scenario | Throat pressure | Maximum along fillet | Max−throat |
|---|---:|---:|---:|
| 60° | −1.080 kPa | +0.056 kPa | 1.136 kPa |
| 95° | +0.188 kPa | +0.895 kPa | 0.706 kPa |
| 120° | +1.080 kPa | +1.463 kPa | 0.383 kPa |

For 60°, the *plan-view* meniscus curvature changes sign at
\(\phi=30°\), i.e. \(x=R\sin30°=10\) µm into the nominal 20 µm
fillet, or arclength \(\ell=10.47\) µm. It does NOT imply the
total 3D capillary pressure changes sign there.

Note that the model predicts a **pressure RISE along the fillet**,
whereas the experimentally extracted quantity is a **pressure DROP
after a burst**. These are different observables, so numerical
similarity of 0.383 kPa (predicted rise at 120°) and 0.386 kPa
(measured candidate drop) is **not quantitative validation**.

## External validation blockers and protocol

1. Map video frames to pressure logger timestamps by **actual metadata
   or a jointly observed event**. TIFF filenames contain original
   camera frame index ranges, but camera and sensor time origins are
   unverified. The metadata-only index script does not download
   any large TIFF and cannot itself prove synchronization.
2. Recover actual width-versus-x coordinates from SEM/video, including
   curvature and roughness, and the independently measured advancing
   contact angles for all four walls from experimental SI.
3. Solve the one-front 3D Young–Laplace free-boundary problem in
   *that measured geometry*, including spatially heterogeneous wall
   shapes. Evaluate stationarity and stability versus wetted volume.
4. Predict independently measurable peak pressure *before* burst and
   at least two post-burst interface locations, with a pressure sensor
   offset, viscous flow and compliance treated as independently
   determined nuisance parameters.
5. Compare with a rectangular constant-angle baseline and the authors'
   own simulations. Publish all signed errors and failed cases.

**Professor interpretation:** the model adaptation and retrieval
protocol are substantive steps, but not a Fuentes-style validation.
The 2025 paper itself already explains curvature reversal, inlet
compliance, and wetting-induced barriers. **Journal readiness 5/10
(unchanged)**. A comparison without independently measured wall shapes
and advancing angles would inflate credibility rather than add it.

## Run (GitHub Actions, manual or relevant PR only)

\`python cah-statics/empirical/audit_darus.py\`
\`python cah-statics/empirical/yan_entry_forward.py\`
\`python cah-statics/empirical/yan_video_index.py\`

Only small reports, CSV and SVG are archived; no raw TIFF.


## Discovered Q1 video filenames: frame-index ≠ logger time

Archive metadata-only inspection recovers:
- 60°: \`60_Q1_water_air_3330-3710_TD.tif\` (frames 3330–3710).
- 95° in the article: **\`90_Q1_water_air_15006-19853_TD.tif\`**.
  This nominal **90 versus 95** discrepancy needs clarification;
  the archive filename must not automatically be renamed to 95°.
- 120°: \`120_Q1_water_air_4700-5357_TD.tif\`.

At 30 frames/s, the middle video's absolute frame indexes imply
500–662 seconds from *the camera's* presumed zero, yet the 95° pressure
recording lasts only 285 seconds. Therefore **camera frame numbering
and pressure logger elapsed seconds do NOT share a proven time origin**.
An offset/synchronization event is essential. The 120° video window
156.7–178.6 s happens to include the 168 s pressure candidate, but
this numerical correspondence is NOT proof of the shared clock.

\`yan_tiff_range_probe.py\` tries to determine whether the 329MB
120° video can be accessed via safe HTTP byte ranges, without
downloading any large image. If range access is unavailable, stop
and request/source a tiny authorized clip for real image-based alignment.


## New: quantitative cross-modal microscope-curvature test (18 October 2026)

A bounded sparse TIFF extractor downloaded original 16-bit image frames
from the **120°** Q1 experiment at only ~10.7 MB of range requests from
the 345 MB original video; images remain attributed to Yan et al.
[DaRUS CC BY 4.0](https://doi.org/10.18419/DARUS-4683).
The camera frames remain precisely labelled with *original* indices.

A central pore darkness proxy from the native images rises sharply:
frame 4900 (index 200) ~0.118; 4905 ~0.543; 4910 ~0.955;
4920 ~0.999. These optical numbers come from an illustrative
central 175×160-pixel region and a fixed raw 16-bit threshold.
They are **not independently calibrated pore-volume saturation**.

### An independent pressure–curvature reconstruction

An OpenCV Canny edge mask and robust least-squares circular fits to the
plan-view menisci (no pressure in image fitting) give:

| Independent image measurement | Value |
|---|---:|
| Original video frame 4900 (index 200) | radius **120.487 µm** |
| Original video frame 4920 (index 220) | radius **429.732 µm** |
| Water–air tension assumed | 0.072 N/m |
| Curvature-derived signed pressure change magnitude | **0.4300 kPa** |
| From candidate measured pressure event (166.15→169.13 s on sensor) | **0.3863 kPa** |
| Central difference | **0.0438 kPa** (11.3% of signal decrease) |
| Fit sensitivity across 9 edge/initial-band choices | **0.394–0.435 kPa** |
| Quoted original pressure-sensor reproducibility | **±0.2 kPa** |

The classic *planar curvature difference* follows from

\[
|\Delta P_{\rm cap}|\simeq
\gamma\left(\frac{1}{R_{4900}}-\frac{1}{R_{4920}}\right),
\]

**only if the out-of-plane curvature and flow/storage corrections
change negligibly between the two images**. This cancellation of
unknown depth curvature can make differences more identifiable
than absolute pressure thresholds.

**CRITICAL: THIS IS NOT A BLIND PRESSURE PREDICTION.**
The radii are **measured from the same experiment**, rather than
predicted from lithographic topography and Young angle; the
classical Young–Laplace relation is well-known. The two image frames
span ~0.667 seconds while the selected sensor release spans ~3
seconds, and camera-logger time-origin offset is **not independently
known**. These cannot yet be called the two endpoints of the *same*
measured event. The numerical agreement could be coincidental.
The sensor ±0.2 kPa precision also precludes a strong fit statistic
from this single provisional pairing.

The source TIFF includes clear films at frames 4900, 4905, 4910,
and 4920; the immersion front traverses the pore rapidly. The
pressure trace's biggest candidate decrease is **later on a naive
shared-zero time axis**, confirming clock synchronization is needed.

**Publication interpretation:** A reproducible cross-modal *candidate*
validation of classical Young–Laplace, a major practical improvement
over graph-digitization. It is neither novelty evidence for our
isospectral two-front pitchfork nor a demonstrated quantitative
prediction of jump positions. **Readiness remains 5/10** until clock
synchronization, measured advancing wettability, actual wall profile,
and a parameter-free forward model produce multiple event predictions.

Reproduce after pressure parsing and sparse image extraction:

```bash
python cah-statics/empirical/yan_crossmodal.py
```

The script generates `yan_crossmodal_provisional.json`,
`yan_radius_fit_robustness.csv`, and an explanatory PNG/SVG figure.
No external video, experimental pressure trace, capillary contact angle,
or pressure sensor offset is fitted to improve the result.
