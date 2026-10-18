# Round 12: sparse 95° original-image recovery and a pressure-clock provenance warning

**Status (18 October 2026): independently observed fast optical front; absolute 3D volume and pressure–volume pairing NOT validated. Journal-readiness 5/10, unchanged.**

## Recovered experiment without a bulk 2.54 GB download

[Yan et al. (2025)](https://doi.org/10.1016/j.jcis.2025.137884), [original archive](https://doi.org/10.18419/DARUS-4683) (CC BY 4.0), DaRUS file ID **351440**, archive MD5 **cc6afebdb4acba25cb94c0f71d15f97b**, **2,544,624,129 bytes**. The **archive filename** is `90_Q1_water_air_15006-19853_TD.tif`, but the actual embedded ImageJ/Nikon metadata contains:

- `Location = /media/lifei1/confocal_experiments/signle_pore_experiments/sorted/water_air/Theta_95/Q1/dranaige/no2/test1035.nd2`.
- `Series 0 Name = test1035.nd2 (series 1)`.
- `dCalibration = 2.485922277608956` (µm per pixel, agreeing with the published ~2.4859).
- `images=4848`; `finterval=0.033328503370285034` s; frame labels explicitly run from `t:15006/21292` to `t:19853/21292`.

This **resolves the nominal 90/95 naming discrepancy for the physical image source**, while not proving that the adjacent 95° pressure file comes from the same experimental run.

A strict HTTP **206 Content-Range** audit established 16-bit Big-Endian 512×512 ImageJ TIFF storage, an uncompressed contiguous image region with exactly 4,848 frames and the next IFD immediately following it. The [source script](sparse_third_wettability.py) reads only header, ~1.63 MB source ND2 info, frame 0 and 48 late frames (**approximately 27.3 MB**, not 2.54 GB) and source pressure header. The original whole-file MD5 is *listed from the archive, not independently recomputed from sparse reads*. No source video was copied into GitHub, and no GitHub Actions workflow was dispatched.

## Reconstructed burst: actual 2D measurement versus conditional liquid volume

A paired-wall geometric detector found the pore box **(116,69)–(438,391)** in pixels at the predeclared nominal 322-pixel width. Five intensity-contrast masks and independent center-band column-wise front arrival were computed. The reference mask is 180 image-intensity counts / 50% per-pixel local dry-to-final contrast; a static exterior is used to correct common illumination drift. A three-pixel wall strip is excluded. This reuses the pre-existing 60°/120° segmentation methodology without calibrating a new threshold to achieve a preferred result.

| Predeclared mask | Start | Final | 10–90% duration, relative to observed mask change |
|---|---:|---:|---:|
| 100-count / 40% | 0.04399 | 0.89351 | 0.20891 s |
| 180-count / 40% | 0.03783 | 0.86740 | 0.21006 s |
| **180-count / 50% reference** | **0.03499** | **0.86740** | **0.19003 s** |
| 180-count / 60% | 0.03282 | 0.86740 | 0.16878 s |
| 260-count / 60% | 0.02809 | 0.84268 | 0.16944 s |

The reference mask crosses 10% and 90% of *its own observed dynamic change* at clip times **160.596679 s** and **160.786704 s**. These are **not** 10% and 90% of independently measured 3D pore liquid volume, unlike a physically calibrated saturation curve. Relative threshold sensitivity is **0.169–0.210 s**, and the temporal sampling itself is **0.03333 s**. The last-frame reference mask recognizes only 86.74% of the nominal projected pore despite visually identifiable corner gas; thin films, near-wall artifacts and 3D height remain unresolved.

The center-band interface front (majority of eligible rows darkened) moved from `x=333` at source clip frame 4822 to `x=395` at frame 4824. This is **62 pixels ≈154.13 µm in 0.066657 s**, or **2.312 mm/s**. The 40%, 50%, and 60% center-band thresholds give the *same* positions for this interval. This is observed **planar front speed**, not a 3D volumetric flow measurement. Relative to a nominal Q1 pump superficial speed of **0.0289 mm/s** over an 800×120 µm pore cross-section, the peak two-frame front speed is about **80-fold faster**. Its physical volume implication still requires a measured transverse/depth-filling field.

[The lightweight measured even-frame table](derived/third_wettability_95_even_frame_audit.csv) and [comparison figure](figures/third_wettability_95_burst.svg) overlay the optical phase mask and normalized front position **without pressure time shifting**. The [standalone stdlib SVG generator](plot_third_wettability.py) reproduces the figure from the frozen table.

**Conditional mass balance only:** Assume every newly dark pixel represents fully filled 120 µm depth. Then, over the reference mask 10–90% duration, 0.8 × (0.86740 − 0.03499) × 0.0768 µL = **0.05114 µL** liquid-equivalent increase, versus nominal pump delivery **0.000528 µL**. The net **~0.0506 µL** excess is of the same order as earlier 60°/120° calculations, but *this numerical resemblance is not independent evidence of common compliance*: comparing fixed fractions of the same nominal pore volume mechanically encourages similar volume differences, and the third case lacks a calibrated 3D saturation.

## Why pressure pairing must be withheld

The Nikon info record contains `dTimeAbsolute = 2459529.1673781597`. **If** this denotes the original ND2 acquisition start as a standard Julian date, conversion yields **2021-11-10 16:01:01.473 UTC**. The matching archive-label pressure file `95_Q1_pressure_water_air.txt` (DaRUS **351431**) begins with sensor `[Context] Date = 11/24/2021 4:34:24 PM`, i.e. **14 calendar days later**, regardless of ordinary timezone differences. The camera file is a subset of a longer 21,292-frame source; the global source frame labels do not remove this interdate discrepancy.

Under that conventional interpretation, these particular *deposited objects* cannot have been recorded on a common hardware clock in the same continuous run. **But `dTimeAbsolute` has not been independently confirmed against the full original ND2 file or a microscope run log**, and might reflect inherited metadata. Treat this as a serious **conditional provenance warning**, not a definite accusation of mislabeling or a contradiction of the published authors' synchronization method.

An **independent** examination of the 95° pressure text shows 17,835 valid samples over 285.36 s at 62.5 Hz. With a 1 s Gaussian smoothing kernel, its strongest negative 1 s pressure change is about **−0.139 kPa** at **257.38 s** after sensor-file start, whereas the video optical burst is near **160.7 s** after *cropped TIFF* start. **No shift was fitted**, because the clock origins cannot be equated. A `+96.7`-second overlay shift chosen from these extrema would be unjustified.

The earlier 60° and 120° TIFFs preserve ImageJ slice labels but **not** the detailed ND2 `info` block. Their deposited pressure texts begin **16 November 2021** and **15 November 2021**, respectively. Thus we cannot run this absolute-date cross-check for them from the cropped TIFFs alone.

## Scientific interpretation and limits

- **Actual new evidence:** independently recovered and time-resolved a third pore-filling front; resolved an archive filename mismatch using original acquisition-folder provenance; identified a concrete potential camera–pressure non-pairing. The fastest observed front is an experimentally based velocity, not a model prediction.
- **What this is not:** a discovery of Haines bursts (already known), a demonstration of new wettability physics, a pressure-extremum fit, a 3D liquid volume measurement, proof of universal ~0.05 µL compliance release, or a falsification of the authors' level-set simulations.
- **Prior art:** Yan et al. already image-track and model pore entry, comparing against established rectangular-channel/Crisp–Thorpe curvature approximations. Fold/pitchfork snapping and capillary instabilities have precedents in [Wells et al. (2018)](https://doi.org/10.1038/s41467-018-03840-6) and [Bisquert et al. (2024)](https://doi.org/10.1021/acs.langmuir.4c02908). The additive two-front energy identity in [our separate theory note](../theory/PRESSURE_VOLUME_TRANSFER.md) is not directly a one-front 95° pressure prediction.
- **Code validation:** [five network-free tests](tests/test_sparse_third_wettability.py) passed, covering ImageJ IFD contiguity, 95° metadata/date parsing, paired-wall geometry, event interpolation and a synthetic moving front. Full source-script replay of the original sparse byte ranges passed, reproducing the pore box, **0.1900254002 s** reference duration, **2.3122220 mm/s** speed, `Theta_95/Q1` source path and pressure header date. SVG parsed as XML and contained both measured traces.
- **Compute discipline:** no full TIFF download, no source-video redistribution, no workflow dispatch; development remains in draft PR #14 on a separate branch.

**Professor's judgment:** mathematical rigor **8/10** (existing conditional identity, not altered); novelty **4/10** (provenance and kinetics are useful but expected capillarity); experimental validation **4/10** (third raw image sequence added, neither 3D volume nor pressure correspondence validated); potential impact **6/10**; submission readiness **5/10 unchanged**.

**Strongest reviewer objection:** The 95° image/mask is another Haines burst already described by Yan et al.; without a trustworthy matching pressure clock and independently determined top/bottom wettability or 3D geometry, no novel predictive Young–Laplace comparison is possible.

**Decisive next action:** Obtain the original `test1035.nd2` time log and a jointly timestamped pressure trace for the **same 95° run**, or a new hardware-triggered synchronized camera/pressure experiment including depth-resolved phase volume and measured fillet. Pre-register intrinsic advancing/receding angles, corner geometry and a single-front 3D Young–Laplace pressure-extremum prediction before looking at the pressure trace. **Journal path:** *Langmuir* or *Physics of Fluids* after that gate; *Physical Review Fluids* requires a genuinely predictive general principle.
