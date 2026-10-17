# Independent experimental benchmark candidates (17 October 2026)

Goal: find observations which distinguish **our 3D two-front statics**
from existing pinning models rather than merely showing stick-slip.

| Source | Measurements | Geometry fit | Current evidence |
|---|---|---|---|
| Poulose et al. 2023 (arXiv:2210.16962; IEEE) | Two interface snapshots, normalized outer L and inner l; ~1mm × 0.38mm; 300μm corrugations, ~10% width | **Direct two-front** geometry; reported plate angles ~90° | Qualitatively consistent with piecewise contact-line steps, smoother inner separation. No independent wall trace, contact-line trajectories or Young sidewall angle. Cannot score jump MAE. |
| Yan, Müller et al. 2025 (JCIS doi:10.1016/j.jcis.2025.137884) | 60/95/120° wettability, single throat/pore front, measured pressure and TIFF images | **Different one-front geometry**; needs a separate finite-depth boundary-value problem and apparatus compliance | Three pressure traces in CC-BY DaRUS (doi:10.18419/DARUS-4683); archive audit script tests access to real raw data. |
| *Towards plug and play filling...* (Biomicrofluidics 2014 PMC4192738) | 3 burst pressures 348±23 / 538±11 / 771±11 Pa | Sharp-edge valves, not sinusoid; single interface | Authors' pre-existing analytical 279/544/809 Pa. Real quantitative standard, **NOT predictions of our theory**. |
| *Microfluidic chain reaction...* (Nature 2022 doi:10.1038/s41586-022-04683-4) | Measured valve thresholds in ~10Pa pressure steps, imaging | Sharp valve geometry; single front | Promising independent valve benchmark; needs independent dimensions/angles. |
| Wells et al. 2018 (Nat Commun 9:1380) | Reproducible droplet snaps on wavelength 2mm, amplitude 0.2mm smooth topography; videos | **Sessile droplet over shaped base**, not confined two-front microchannel | Closest mechanistic prior art; cannot directly plug into our pressure-volume solver. Authors already model snap bifurcations. |
| Al Balushi et al. 2023 (Langmuir doi:10.1021/acs.langmuir.3c00914) | Several mixed stick-slip regimes on pillars, water/ethanol | Chemistry/mixture and pillared free droplets | Weak direct applicability. |

**Critical boundary:** none of these independently verifies the new 3D
isospectral phase-sensitive bifurcation sign result. The Poulose experiment
has sidewall waviness but does NOT provide a manufactured pair with identical
Fourier magnitudes and opposite phases. Its apparent 30° cannot be used
as an independently measured Young sidewall angle. Yan is one moving
front entering a throat; it cannot validate the two-front pitchfork.

## Real-data audit

Run the manual or PR-scoped **CAH experiment — DaRUS pressure audit**
Action. It queries the original Utrecht/Stuttgart dataset metadata and
three lowest-flow water–air pressure files, and emits *only a small
derived JSON structure/quality report* (file checksums, column format,
sample header, numeric ranges). No large TIFF download, raw data
redistribution, ongoing schedule, or assertion of experimental fit.

The source paper's 2025 low-flow behaviour is a **positive control for
geometric curvature reversal** at 60° but not at 95° or 120°; the
experimenters already describe this mechanism. Prediction of *pressure
thresholds* requires independently measured pore-throat profile,
four-wall contact angles, apparatus pressure calibration and compliance.
The old rough figure readings do NOT count as raw experimental data.

## The next genuine model test

1. Independent optical metrology of a smooth corrugated channel,
   contact angles on all 4 walls, and each meniscus' position over time.
2. Predeclare predicted pressure extrema, volumes and their uncertainty
   intervals for the measured geometry; leave one experiment blinded.
3. Compare baseline models: (a) uniform rectangular capillarity, (b)
   conventional 2D geometric model, (c) current 3D capillary-energy model.
4. Report MAE for jump positions/volumes/pressure and failure cases;
   fit no wall profiles or contact angles to observed jumps.
5. The phase-only topography result needs a deliberately fabricated
   isospectral pair or a measured equivalent; other experiments are
   **general mechanism checks**, not validation of that new claim.
