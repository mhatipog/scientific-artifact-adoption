# Independent Yan (2025) pore-volume audit

The source paper and raw experimental TIFFs are [Yan et al., JCIS 2025](https://doi.org/10.1016/j.jcis.2025.137884) and [DaRUS-4683](https://doi.org/10.18419/DARUS-4683), CC BY 4.0. No raw microscope videos are committed.

- [Full measured results and scientific limitations](ROUND11_PHASE_VOLUME.md)
- [Vector figure comparing published saturation and our independent masks](figures/phase_saturation_vs_published.svg)
- Derived lightweight trajectories: derived/segmented_saturation_10hz.csv and frozen 30-Hz crossing summaries in derived/measurements.json.

Reproduce **offline** with Python 3.11+ and numpy, scipy, Pillow:

    python tests/test_phase_mask.py
    python plot_comparison.py

Reprocess one original experiment **only when deliberately requested** (200-345 MB, source MD5 verified):

    python segment_original_tiffs.py --angle 60 --download-original
    python segment_original_tiffs.py --angle 120 --download-original

No scheduled workflows are added or triggered. This subfolder is intentionally outside legacy experiment-path PR workflows that would otherwise re-download pressure files and consume GitHub Actions quota.

The geometry detector is specialized to these exact two TIFF stacks. A two-dimensional water mask times nominal 120-µm depth is a hypothesis, not a measured 3D liquid volume or capillary pressure.
