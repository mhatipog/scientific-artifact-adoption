# Independent Yan (2025) pore-volume audit

The source paper and raw experimental TIFFs are [Yan et al., JCIS 2025](https://doi.org/10.1016/j.jcis.2025.137884) and [DaRUS-4683](https://doi.org/10.18419/DARUS-4683), CC BY 4.0. No raw microscope videos are committed.

- [Full measured results and scientific limitations](ROUND11_PHASE_VOLUME.md)
- [Vector figure comparing published saturation and our independent masks](figures/phase_saturation_vs_published.svg)
- Derived lightweight trajectories: derived/segmented_saturation_10hz.csv and frozen 30-Hz crossing summaries in derived/measurements.json.

Reproduce **offline** with Python 3.11+:

    python -m pip install -r requirements.txt
    python tests/test_phase_mask.py
    python plot_comparison.py

Reprocess one original experiment **only when deliberately requested** (200-345 MB, source MD5 verified):

    python segment_original_tiffs.py --angle 60 --download-original
    python segment_original_tiffs.py --angle 120 --download-original

No scheduled workflows are added or triggered. This subfolder is intentionally outside legacy experiment-path PR workflows that would otherwise re-download pressure files and consume GitHub Actions quota.

The geometry detector is specialized to these exact two TIFF stacks. A two-dimensional water mask times nominal 120-µm depth is a hypothesis, not a measured 3D liquid volume or capillary pressure.


## Sparse third-wettability / pressure-clock audit (Round 12)

The archive labels its 2.54-GB image `90_Q1_water_air`, but embedded Nikon metadata explicitly identifies the original `Theta_95/Q1` folder. We used HTTP 206 byte ranges to read only ~27 MB, confirm the contiguous 4,848-frame TIFF, recover a 0.169–0.210 s optical burst and a front speed of ~2.31 mm/s. **The new 95° pressure file has a potentially incompatible absolute acquisition date**; no camera–pressure shift was fitted.

- [Independent measurement, metadata and clock interpretation](ROUND12_SPARSE_95_CLOCK.md)
- [Sparse original-TIFF + provenance executable](sparse_third_wettability.py)
- [Offline synthetic/header/time/geometry tests](tests/test_sparse_third_wettability.py)
- [Frozen even-frame raw-image-derived event samples](derived/third_wettability_95_even_frame_audit.csv)
- [Optical mask / physical interface comparison](figures/third_wettability_95_burst.svg) and [reproducer](plot_third_wettability.py)

Reproduce offline without any downloads:

    python tests/test_sparse_third_wettability.py
    python plot_third_wettability.py

Opt into a *manual* approximately 27-MB source byte-range audit (not an Actions job):

    python sparse_third_wettability.py --download-ranges

Source TIFFs and proprietary parent `.nd2` files are not committed. Reading byte ranges cannot independently verify the supplied full-file MD5. **No 3D volume calibration or 95° image/pressure hardware synchronization is claimed.**


## Loss-guarded pump-only depth certificate (Round 13)

The [Round-13 proof and original-image results](ROUND13_CRITICAL_DEPTH.md) invert mass conservation to ask what **independently measured incremental water thickness** would falsify pump-only filling without requiring a camera–pressure time shift. In the clean 95° 4822→4824 interval, five frozen optical masks yield **0.983–1.400 µm** worst-case critical depth, even if all reversed-class pixels release the full nominal 120-µm pore depth. The earlier window is **not** segmentation-robust. This is a *conditional falsification threshold*, not a measurement of 3D liquid volume, reservoir release or pressure.

- [Sparse source-image classifier and activation/loss balance](critical_depth_certificate.py)
- [Ten frozen window/threshold measurements](derived/critical_depth_95_windows.csv)
- [Reproducible SVG comparison](figures/critical_depth_95_guarded.svg) and [offline plot generator](plot_critical_depth.py)
- [Seven zero-network regression tests](tests/test_critical_depth_certificate.py)

Reproduce the figure and tests **without network**:

    python tests/test_critical_depth_certificate.py
    python plot_critical_depth.py

Opt in only when deliberately re-reading ~4.72 MB of the archived 2.54-GB TIFF:

    python critical_depth_certificate.py --download-ranges

**The still-missing decisive input is calibrated *incremental* liquid-column thickness, not optical intensity, and pressure synchronization remains unverified.** This research PR remains draft and publication readiness remains 5/10.

## Retained-wet optical response audit (Round 15)

The [Round-15 original-frame diagnostic](ROUND15_RETAINED_OPTICS.md) shows that even after excluding a 7-pixel band around the classified wet/dry boundary, the original 95° frames 4822→4824 contain **33,175 persistently wet pixels** with median normalized optical change **+0.16759**, compared with **+0.00676** in stable-off pixels. Five frozen classifiers and 15-pixel boundary exclusion preserve the qualitative observation. These are **optical values, not calibrated water-column depths or volumes**; they expose a limitation of the binary area measure without resolving the mass budget.

- [Manual strict-206 original TIFF audit](optical_retention_audit.py) (four source frames only; emits 20 segmentation/erosion rows)
- [Five-row 7-pixel-boundary derived CSV](derived/retained_optical_95_4822_4824_7px.csv)
- [Four offline tests](tests/test_optical_retention_audit.py)

Run tests offline:

```sh
python cah-statics/experimental_volume/tests/test_optical_retention_audit.py
```

To independently repeat the **manually opted-in** public-TIFF read (~2 MiB, original media never committed):

```sh
python cah-statics/experimental_volume/optical_retention_audit.py --download-ranges
```

The calibration of **incremental 3D water thickness** and camera–pressure synchronization remains missing. Both experimental and theoretical PRs stay draft; publication readiness stays **5/10**.
