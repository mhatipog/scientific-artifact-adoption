# Round 10: Conditional constant-compliance holdout (18 October 2026)

**Status: an experimental consistency check, NOT verification of our capillary-statics theory.**

Yan et al. (JCIS 2025, doi:10.1016/j.jcis.2025.137884) pressure and video data are available under CC BY 4.0 at doi:10.18419/DARUS-4683. We registered each original sensor waveform to published Figure 5, and each original TIFF optical filling signal to published Figure 4. **These are two separate time-axis fits to experimental plots, not a shared hardware trigger.** No pressure amplitude, capillary parameter or Young angle is fitted.

| Conditional 10–90% optical event | 60° | 120° |
|---|---:|---:|
| Optical duration | 3.488 s | 1.114 s |
| Pump delivery | 0.00969 µL | 0.00309 µL |
| Implied liquid volume gain (assumed optical-to-volume calibration) | 0.06144 µL | 0.05652 µL |
| Estimated net storage release | 0.05175 µL | 0.05343 µL |
| Pressure change on registered experimental time axis | −0.2720 kPa | −0.0753 kPa |
| Inferred effective compliance | 0.1902 µL/kPa | 0.7099 µL/kPa |

**Held-out test:** calibrating a single constant compliance from the 60° experiment predicts a −0.2809 kPa pressure change for the 120° experiment, versus −0.0753 kPa conditionally inferred from its recorded pressure waveform. The prediction-to-observation magnitude ratio is 3.73. This is an internally computed disagreement of a simple apparatus **null model**, not a claim that the capillary theory is inaccurate or accurate.

Additional camera time shifts of ±0.4 s yield apparent compliance ranges **0.190–0.199 µL/kPa** (60°) and **0.359–1.103 µL/kPa** (120°). These are **sensitivity scenarios, NOT confidence intervals**. The ranges do not overlap *under that assumed time-window test*.

The inference is **not physically identifiable yet**. Optical intensity is not independently calibrated as pore volume; the 120° experiment has residual corner gas; pressure recorded at the inlet can also contain hydraulic resistance and compliance dynamics; and the synchronization relies on digitizing published plots. Independently assuming 50%–100% of the nominal optical-volume conversion for each case still does not restore equal compliance, but that conversion range is itself hypothetical, not a measured bound. Do not report the difference as evidence that apparatus compliance is nonlinear.

**What changes research direction:** Segment the pore liquid area/volume in each raw TIFF frame, including corners, calibrate depth and advancing angle, and obtain genuine camera–sensor acquisition metadata. After that, independently predict the pressure extrema with the *single-front* 3D Young–Laplace solver. Yan is not a direct test of our distinctive two-front isospectral pitchfork.

**Professor assessment — this round:** experimental methodology improved; direct independent capillary-pressure prediction still absent. Journal submission readiness **5/10, unchanged**. Prospective journals Langmuir or Physics of Fluids if a model beats prior baselines on multiple held-out events.

To reproduce: run bridge_camera_pressure.py, then holdout_compliance.py using the frozen independent reference CSVs and recorded TIFF-derived optical timing summary. GitHub Actions generates the PNG/SVG and JSON evidence.