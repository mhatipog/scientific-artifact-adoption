# Round 11 (18 October 2026): original TIFF phase masks and volume identifiability

**Status: independent optical segmentation; method-sensitive quasi-2D volume proxy, not a validated capillary-pressure prediction. Submission readiness remains 5/10.**

## Provenance and source geometry

- L. Yan, J. C. Müller, T. van Noorden, B. Weigand and A. Raoof, *J. Colloid Interface Sci.* **696** (2025), 137884, [DOI](https://doi.org/10.1016/j.jcis.2025.137884), CC BY 4.0. Original TIFF data [DaRUS-4683](https://doi.org/10.18419/DARUS-4683), CC BY 4.0. **Original videos are not redistributed.**
- 60_Q1_water_air_3330-3710_TD.tif: 381 frames, MD5 3dbfc99446d61b89c905057cd360d5b7. 120_Q1_water_air_4700-5357_TD.tif: 658 frames, MD5 1beeb40ae47b8c7299aba034e40c5c5f. Both original MD5s independently checked before reading **all 1,039 frames**.
- Nominal pore 800 × 800 × 120 µm; pixel pitch 2.4859 µm, depth assumed uniform; nominal pore volume 0.0768 µL; Q1 pump 0.01 mL/h = 0.00277778 µL/s.
- Geometry conflict remains: Fig. 1d SEM caption calls the rounded pore opening about **20 µm** radius; Sec. 4.2 gives **5.6 µm** for an operational contact-line pinning region. These are not a uniquely specified fillet profile. A 5.6 µm feature is only **2.25 image pixels**, while 20 µm is **8.05 pixels**. Do not infer 3D corner curvature from these videos.

## Direct measurements and reproducibility

The manually executable segment_original_tiffs.py downloads just one original TIFF when explicitly requested, validates its source MD5, locates a pair of wall lines separated by the nominal 322-pixel pore span, and classifies water by initial-to-final local contrast at every image frame after exterior illumination correction. A three-pixel wall strip is deliberately excluded. Five predetermined contrast settings explore method sensitivity. The method does **not** label every non-water pixel as gas: such pixels may include residual corner air, poorly resolved films, dust or wall artifacts.

The detected pore boxes (x_min, y_min, x_max, y_max) are **(96,62,418,384)** for 60° and **(96,31,418,353)** for 120°. The detector's search windows are deliberately specialized for these two 512px datasets. The reference mask uses a 180-intensity-count contrast threshold and 50% of each pixel's reference darkening. Results below are computed at full TIFF frame rate, then the plotted derived curves are thinned to approximately 10 Hz to keep the public repository small.

| Quantity | 60° | 120° |
|:--|--:|--:|
| Frames | 381 | 658 |
| Final mask fraction, permissive (100 counts/40%) | 94.961% | 92.163% |
| Final mask fraction, reference (180 counts/50%) | 93.967% | 91.232% |
| Authors' published final saturation | ~100% | ~92–92.5% |
| Mask first crossing of 10% physical target | 6.850 s | 6.773 s |
| Mask first crossing of 90% physical target | 11.029 s | 8.225 s |
| 10–90% duration | 4.178 s | 1.452 s |
| Conditional pore liquid-volume gain | 0.06144 µL | 0.05652 µL |
| Nominal pump delivery in interval | 0.01161 µL | 0.00403 µL |
| **Net storage release under mass balance** | **0.04983 µL** | **0.05249 µL** |
| Pore filling/pump volume ratio | 5.29 | 14.01 |

Target percentiles are absolute fractions of **nominal pore volume**, 0.1/0.9 at 60° and 0.092/0.828 at 120° using the paper's independently reported final ~92%. No masks were scaled to force initial/final saturation agreement.

**Positive check:** The 120° permissive mask gives 92.16%, within ~0.34 percentage point of the digitized published ~92.5% terminal saturation, without adjusting the mask to those data.

**Negative check:** Even the permissive 60° mask ends at 94.96% although the authors observed a completely filled pore. Therefore about 5% remains unsegmented; physical films/wall-adjacent pixels, static optical artifacts and mask omission have not been independently separated.

## Comparison with the authors' saturation and earlier brightness analysis

The [self-contained comparison figure](figures/phase_saturation_vs_published.svg) overlays the full raw-image mask trajectories and the authors' independently digitized Fig. 4 without optimizing **any** relative time shift.

Against the published Fig. 4 saturation at identical printed clip time, reference-mask RMS errors over t=2–12 s are **5.16 and 4.15 percentage points** (60°, 120°). Early water is systematically missed. At 5 s the publication gives ~8.0% and ~6.5% versus 2.2% and 0.9% from our 60° and 120° reference masks. In the main bulk-invasion interval the agreement is often much closer. These curve errors exceed the authors' nominal ~2% saturation uncertainty. They prevent a publication-grade direct physical-volume claim.

Reanalyzing the authors' own segmented Fig. 4 yields 10–90% durations **4.920/1.865 s** and conditional storage releases **0.04777/0.05134 µL**. Both sources imply ~0.05 µL fluid beyond direct pump delivery, but this is a volume-budget inference, **not** direct bubble/compliance measurement. The authors themselves already identified a possible inlet gas bubble as the mismatch mechanism (Sec. 6).

## Pressure identifiability: corrected null-model assessment

The earlier brightness-proxy constant-compliance ratio C120/C60 ≈ 3.73 does **not** survive a better optical-to-volume conversion. Using authors' Fig. 4 physical saturation and the frozen matched source pressure, and **conditionally assuming** publication Fig. 4 and Fig. 5 have identical experimental time origins, the pressure changes are −0.316 and −0.151 kPa and the apparent compliances **0.151 and 0.339 µL/kPa**, ratio **2.25**.

A shared C calibrated on 60° predicts approximately **−0.340 kPa** for 120°, versus the conditionally observed −0.151 kPa, a ~**0.188 kPa** difference. Yan et al. report sensor system accuracy including reproducibility **±0.200 kPa**, not a specified error bound on the difference of readings in one run. Differential bias of the needed scale is not independently excluded; conversely a common offset cancels, so the specification does not prove such a bias. The common time-origin hypothesis also lacks a verified shared hardware trigger in the deposited TIFF metadata. Accordingly **the shared constant-C hypothesis is neither falsified nor confirmed**.

A perfect planar mask still cannot measure true 3D pore water volume without the depth-dependent filling height h_w(x,y,t): V = integral over the pore of h_w dA. Taking h_w=120 µm is an assumption. This and the underresolved 5.6 µm corner block a blind 3D pressure-extremum prediction from these videos alone.

## Relative novelty and rigor

The supplementary authors' level-set and Crisp–Thorpe baseline already achieves ~1.8% mean curvature deviation in the 95° reference case. The present analysis does **not** outperform it. Wells et al. (2018) and Bisquert et al. (2024) already document capillary fold/pitchfork snapping; storage-driven Haines jumps are also known. The separate exact two-front transfer identity is explained in theory/PRESSURE_VOLUME_TRANSFER.md and tested, but is a consequence of an additive symmetric energy, not currently an established new physical mechanism.

**Professor's scores:** mathematical rigor **8/10** (exact conditional theorem); novelty **4/10**; experimental validation **4/10** (raw-mask endpoint check but incomplete films); potential impact **6/10**; current journal-readiness **5/10**, unchanged.

**Strongest reviewer objection:** This is still a method-sensitive uniform-depth optical surrogate, substantially reproducing the authors' already-published saturation, with neither validated pressure correspondence nor a novel predictive result.

**Decisive next measurement:** Independent SEM/profilometer inlet fillet, actual advancing angles on side and plates, a hardware-common-trigger pressure/camera experiment, ideally with z-resolved fluorescence. Then pre-register and solve one-front 3D Young–Laplace, compare *preburst* pressure extrema and transition positions to the authors' strong level-set/Crisp–Thorpe baseline.

**Possible journals after that gate:** Langmuir or Physics of Fluids. Physical Review Fluids would need a demonstrated general predictive principle.

Sources: [Yan et al. 2025](https://doi.org/10.1016/j.jcis.2025.137884), [data](https://doi.org/10.18419/DARUS-4683), [Wells et al. 2018](https://doi.org/10.1038/s41467-018-03840-6), [Bisquert et al. 2024](https://doi.org/10.1021/acs.langmuir.4c02908).
