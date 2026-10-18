# Round 15 — retained-wet optical dynamics are not a 3D volume measurement

**Scope/status:** independent original-image diagnostic; observed bright-field intensity, **not** calibrated liquid-column height, physical water volume, or synchronized pressure. The source is Yan et al., *Journal of Colloid and Interface Science* (2025), doi:10.1016/j.jcis.2025.137884, original deposited acquisition doi:10.18419/DARUS-4683 (CC BY 4.0). The archive filename reads 90° Q1, but the embedded Nikon location identifies `Theta_95/Q1`. No source TIFF is stored in GitHub.

## The falsifiable question

Round 13 counted newly positive and disappearing optical-class pixels to derive an **explicitly conditional** critical incremental water depth of 0.983–1.400 µm for frames 4822→4824. Does the classification conceal additional *optical* evolution in pixels already positive at both frames? A positive result would preclude treating binary-class saturation as a complete description of the optical event; it would **not** quantify unseen water redistribution.

## Original data, frozen algorithm, controls

Read only the TIFF header, source frame 0, final frame 4847, and frames 4822 and 4824 of DaRUS file ID **351440** through strictly verified HTTP 206 ranges (four frames ~2 MiB plus header, no full-video redistribution). The frame interval is **0.033328503370285034 s**; the selected interval lasts **0.0666570067 s**. The original pore-wall detector returns `(116,69,438,391)`, corresponding to the independent Round-12 geometry audit.

Use the previous Gaussian blur `sigma=1`, pore interior excluding three wall pixels, median exterior-field illumination correction, and the identical five frozen pairs of contrast-count and fractional mask thresholds. For pixel `i`, define its **dimensionless optical progress**

\[
r_i(t)=\frac{I_i(0)-I_i(t)+\delta(t)}
            {I_i(0)-I_i(4847)},\qquad
\delta(t)=\operatorname{median}_{\mathrm{exterior}}[I_j(t)-I_j(0)].
\]

Only pixels with strictly positive reference-to-final optical contrast exceeding their predeclared count threshold are eligible. Define `retained-wet` to mean above that threshold in **both** frame 4822 and frame 4824. Define `stable-off` to mean below it in both frames. To suppress optical mixing from a displaced meniscus boundary, retain only pixels whose Euclidean distance from **any not-retained-wet pixel** exceeds **7 pixels** (~17.4 µm). This is an interior-mask erosion, not a confidence interval or physical-depth measurement. A 15-pixel erosion gives the same qualitative conclusion.

## Measured optical response, 7-pixel retained-mask boundary exclusion

| Contrast counts, fraction | Retained-wet interior pixels | Median \(\Delta r\) | Retained with \(\Delta r>0.1\) | Stable-off median \(\Delta r\) | Newly positive pixels | Reversed pixels |
|---|---:|---:|---:|---:|---:|---:|
| 100, 0.4 | 39,767 | 0.185 | 81.2% | 0.0069 | 24,580 | 37 |
| 180, 0.4 | 39,112 | 0.187 | 82.1% | 0.0065 | 23,686 | 22 |
| 180, 0.5 | 33,175 | 0.168 | 80.0% | 0.0068 | 28,528 | 29 |
| 180, 0.6 | 28,063 | 0.151 | 76.8% | 0.0071 | 33,773 | 38 |
| 260, 0.6 | 27,553 | 0.152 | 77.5% | 0.0067 | 32,930 | 20 |

The reference choice (180 counts, 0.5 fraction) has **33,175 retained-wet interior pixels**. Their median normalized optical change is **+0.16759**, and **80.04%** increase by more than 0.1, against a median **+0.00676** in 11,520 stable-off controls. Across five original threshold choices, the retained interior median is **+0.1508 to +0.1868**. After the more severe 15-pixel erosion, it remains **+0.1409 to +0.1815**.

The measured exterior drift in original gray counts changes from about **−2.816** (frame 4822) to **−3.509** (frame 4824); the estimator subtracts its effect. Artificial whole-frame intensity offsets leave the offline regression unchanged. These are deterministic descriptive statistics on a single spatially correlated movie event; the five segmentations are **not** five experimental replicates. No p-value or confidence interval for independent specimens is implied.

## Conclusions that *do* and *do not* follow

**Observed:** Strong, spatially localized normalized optical change persists in the regions classified as wet at both ends of the rapid event. It is not explainable solely by classifying newly crossing pixels or by a thin Gaussian-blurred mask boundary. The previously tabulated binary area alone is therefore an incomplete optical description of the burst.

**Not observed:** We cannot identify whether that retained-pixel signal represents increasing actual water-column depth, meniscus curvature/refraction, refractive-index changes, bubbles, contrast formation, or mixed optical phenomena. In particular, it **cannot** be added to the Round-13 volume budget as measured volume. No pressure trace is synchronized to these frames; the acquisition-date conflict remains unresolved. The 95° nominal sample is not the two-front pitchfork experiment.

**Structural inverse limitation:** Even under the optimistic scalar model `r=f(h)` with an unknown monotone optical transfer function, distinct choices of \(h=f^{-1}(r)\) give the same image contrast but different \(V(t)=\int h\,dA\). For example, after clipping \(r\) to [0,1], \(h_1=d\,r\) and \(h_2=d\,r^2\) obey the same 0 and \(d\) endpoints yet predict different intermediate liquid columns. This construction only demonstrates non-identifiability in the **uncalibrated inverse model**; it is not a claim that either transfer is the microscope's physical law.

Thus the honest next step is an **independent spatially resolved depth calibration** on the actual throat–pore chip, ideally synchronized volumetric/confocal or calibrated interferometric data plus pressure and pump flow on a shared acquisition clock. In particular, test the predeclared Round-13 **area-averaged incremental wet-column thickness threshold of at least ~1.5–2 µm over newly positive pixels**, together with changes over the retained-wet interior. No arbitrary camera–pressure shift is admissible.

## Reproducibility, use, limits

- Executable: `optical_retention_audit.py`, manually invoked with `--download-ranges`; prints the 20 full (five classifiers × four erosions) rows, never downloads automatically or saves original TIFFs.
- Frozen lightweight derived result: `derived/retained_optical_95_4822_4824_7px.csv` (the five classifier rows at 7-pixel exclusion).
- Offline synthetic tests: `tests/test_optical_retention_audit.py` (normalization, control group, boundary exclusion, illumination shift, malformed images).
- No Actions workflow, no pressure registration, no 3D Young–Laplace validation; the experimental PR remains **draft**.
- The observation strengthens the case for depth calibration but **does not raise journal-submission readiness above 5/10**.
