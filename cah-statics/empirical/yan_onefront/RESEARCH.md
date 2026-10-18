# Round 7 — Yan (2025): controlled one-front corner test

## Experimental correction and calculable prediction

The 2025 Yan, Müller et al. paper (JCIS 696, 137884;
https://doi.org/10.1016/j.jcis.2025.137884)
reports a **rounded pore corner radius ~5.6 µm**, not 20 µm.
Its contact-angle interpretation also distinguishes *static advancing angle*
from nominal substrate wettability in the low-flow regimes.

At a locally circular 90° pore corner, the capillary pressure's
plan-view contribution changes as

`P_xy(φ)=2γ cos(θ_adv+φ)/w_local`.

Hence for *assumed* θ_adv=60° the local curvature crosses zero after
φ=30°. The predicted distance **along the actual rounded wall** is

`s_c = R_c (π/6) = 2.93215 µm`.

The corresponding displacement projected along a straight axis is
`R_c sin(π/6) = 2.8000 µm`. In a 95° or 120° advancing-angle
scenario the curvature does not cross zero over the 0°–90° corner
turn. This is a **geometrical consistency test**, not a newly
discovered mechanism: the Yan authors describe wettability-dependent
curvature reversal.

The earlier 10.47 µm result arose from an **assumed** 20 µm corner and
**must not** be presented as a numerical prediction for this experiment.

## Why the pressure signals do not yet validate the theory

Three true low-flow DaRUS pressure records have been retrieved and parsed
at about 62.5 Hz. Instrument output is **kPa**, whereas the internal
calibration entry AlgorithmUnit is **psi**. They contain 56,180,
17,835 and 51,186 valid time-stamped samples at nominal 60°, 95°, 120°.
The strongest five-second negative excursions were ~−0.601,
−0.339 and −0.361 kPa respectively.

Those are **not identified pore-filling pressure jumps**, since the
images have not been aligned to the pressure-sensor clock. They cannot
yet be used as target labels for fitted or blind burst-pressure error
calculations. A capillary pressure threshold requires:
- actual inlet/pore *local wall profile* (not one throat width);
- independently measured static **advancing** angles on relevant walls;
- the proper top-bottom plate curvature / 3D meniscus shape;
- pressure baseline/offset and image-to-pressure time registration;
- consideration of syringe compliance and trapped gas during rapid events.

The accompanying numeric pressure plot assumes γ=0.072 N/m and a
150 µm local width **solely to display the scale of one planar
curvature contribution**. It is *not* a full predicted measured
pressure for any of the three archived traces.

## Findings vs publishing

The geometric prediction is more tightly tied to independently
reported morphology, but **not a publishable new prediction**:
the experiment already recognizes the curvature-reversal physics.
What could be publishable is accurate *blindly predicted* critical
pressure / interface position from measured geometry and wettability,
preferably beyond the paper's own established simulations and the
Crisp–Thorpe approximation the authors validated.

**Submission readiness after this round: 5/10 theoretical program,
2/10 independent experimental prediction, 5/10 overall at most.**
Next: obtain image timestamp synchronization and independent wall
segmentation, and use the existing 3D capillary solver on the
single-front throat/pore, not the two-front isospectral pitchfork.
