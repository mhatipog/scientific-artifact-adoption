# Experimental research — preliminary raw-data comparison, 17 Oct 2026

## Actual acquisition and measured facts

The public DaRUS (doi:10.18419/DARUS-4683; CC BY 4.0)
pressure recordings from **Yan et al. (JCIS 2025)** are now
programmatically accessible in GitHub Actions. The original files are
instrument exports from FISO FPI-HR-2 / FOP-M260, sampled every **0.016 s**
with each data row as date, clock time, pressure, validity.

The instrument metadata explicitly states **\`Unit=kPa\`**. A separate
\`AlgorithmUnit=psi\` entry refers to the internal calibration range,
not the column's exported unit; do not inadvertently multiply the measured
kPa values by 6.895.

| Low-flow water–air experiment | Valid samples | Recording duration | Pressure range (kPa) | Largest 5s decrease (kPa) |
|---|---:|---:|---:|---:|
| 60° | 56,180 | 898.86 s | −2.070 to +0.534 | −0.601 |
| 95° | 17,835 | 285.34 s | −0.411 to +1.056 | −0.339 |
| 120° | 51,186 | 818.96 s | −1.851 to +1.196 | −0.361 |

All three downloads passed comparison of their **original MD5 checksum**
with the Dataverse metadata. The tabulated five-second drops are defined
as the largest decline between consecutive five-second-offset,
one-second-block-averaged pressure values. **They are not uniquely
identified pore-filling pressure jumps.** Inspection of the full curves
shows multiple pressure transitions and extended drifts, especially
for the 60° and 120° series.

![Actual archived measured pressure signals](outputs/empirical/measured_pressure_traces.svg)

The raw TIFF images are **not** downloaded: the larger recordings can
be hundreds of megabytes and event alignment requires matching
individual image timestamps. The pressure-event annotations here are
diagnostic, not the peer-reviewed paper's event definitions.

### Precisely what our CAH statics can explain now

1. **Poulose et al. (2023):** the existing finite-volume two-front
   equilibrium model can create discrete contact-line jumps during
   volume-driven evaporation. Their figure shows stair-stepped outer
   contact-line separation \`L\`, smoother inner-meniscus distance \`l\`,
   approximately 1mm × 0.38mm channel, ~300µm laser corrugation.
   Their top/bottom plate contact angles are reported near **90°**,
   making our neutral-plate benchmark directly relevant.
   **No quantitative agreement yet:** measured wall profile, intrinsic
   sidewall Young angle and individual time-resolved front positions
   remain unavailable. Magnetic changes in evaporation *rate* are
   a different transport question not explained by the static model.
2. **Yan et al. (2025):** the wettability-dependent pore-throat barrier
   and interface-curvature reversal follow geometrically from the
   local wall orientation and Young angle. The authors already present
   this explanation. Our full two-front pitchfork prediction is NOT
   being compared: Yan uses one front and a throat/pore expansion, with
   external syringe pump and possible compliance/gas storage.
   The retrieved pressure data show real pressure transitions in kPa,
   but there is **no parameter-free predicted burst-pressure series**
   from our model for their independently measured geometry.
3. **Hagmeyer et al. (2014), Biomicrofluidics 8:056501:** three
   independent measured stop-valve burst pressures
   348±23, 538±11 and 771±11 Pa, compared in *that paper*
   with **their own** analytical predictions 279, 544, 809 Pa.
   Its published analytical baseline MAE is **37.7 Pa** and mean
   absolute percentage error **8.62%** across these three valves.
   Those are *not predictions from our CAH model*. It is a clean
   quantitative independent benchmark **once** we reconstruct actual
   valve geometry and angle boundary conditions.
4. **Yafia et al. (2022), Nature 605:464:** stop-valve experiments
   with roughly 10Pa applied-pressure increments and finite-element
   simulations; strong modern comparison target for future geometry
   adaptation. It is NOT a two-front droplet test.
5. **Wells et al. (2018), Nat Commun 9:1380:** experiment confirms
   a snap/bifurcation interpretation, but the authors already solved
   it; substrate is a *wavy support for sessile 3D drops*, not our
   two-sided confined channel. It is prior art rather than an
   independent novel explanation.

## Publishability assessment — empirical round

**Theoretical readiness 5/10; empirical validation 2/10; journal-submission
readiness remains 5/10 at most.** Real raw data acquisition and plotted
signals are a meaningful reproducibility advance. However, none of
these retrieved recordings yet verifies the unique contribution
(the *isospectral phase-switch* in freely moving confined menisci).

A journal-worthy next result demands **blind prediction of event thresholds**:
image-aligned pore-pressure extrema for Yan's *single-front* geometry;
or independently measured, time-resolved contact-line positions for
Poulose's *two-front* geometry. The baseline 2014 valves can be a
quantitative control. Do not pool experiments with incompatible
boundary conditions or label a visually similar transition a fitted
match.

All reproducible calculations and unaltered-source citations are in
the repository. The manual/PR-only workflow downloads three pressure
files on an ephemeral runner, exports small summaries/decimated figures,
and does not schedule ongoing downloads or commit raw datasets.


## Round 7: nominal Yan pore entry and raw pressure alignment

See [Round 7 numerical baseline and explicit non-agreement](ROUND7_YAN_FORWARD.md). The ideal 20-µm entrance model is a *reduced pressure closure*, not a full 3D simulation or a fitted experimental match. The manual/PR Actions also index the original 30-fps TIFF frame ranges without downloading them.
