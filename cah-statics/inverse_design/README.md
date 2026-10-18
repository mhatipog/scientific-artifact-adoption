# Guarded fifth-order phase design (Round 14)

This is a **theoretical inverse-design audit** in the existing contact-angle hysteresis programme. It is **not** part of the experimental Yan image-volume reconstruction. The latter remains separately documented in draft research PR #14.

- [Full scientific note: mathematical, numerical and novelty qualifications](ROUND14_GUARDED_PHASE_PAIR.md)
- [Independent analytic-derivative + fold/root/uncertainty executable](phase_pair_window.py)
- [Eight zero-network unit tests](tests/test_phase_pair_window.py)
- [Frozen nominal and stress summary](derived/guarded_pair_summary.json)
- [45-stratum aggregate of 405 geometric/wetting combinations](derived/guarded_pair_45_strata.csv)
- [Deterministic vector figure](figures/phase_pair_robust_window.svg) and [figure generator](plot_phase_pair_window.py)

The selected 5th/7th-harmonic pair retains **exact** Fourier power/extrema/correlations through order five at nominal geometry, while its local reduced-model quartic signs differ. An independent uncertainty grid checks local root isolation, sidewall angle variation, relative Fourier fabrication errors, sampled wall clearance and the **special neutral-plate** case. The design explicitly avoids a nearby pressure-stationary fold creation; numerical robustness over 405 samples is **not interval certification**.

Use NumPy and SciPy. All regular commands are **offline**:

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_phase_pair_window.py' -v
python cah-statics/inverse_design/phase_pair_window.py
python cah-statics/inverse_design/phase_pair_window.py --grid-csv
python cah-statics/inverse_design/plot_phase_pair_window.py
```

No workflow is installed or dispatched. No experimental media are included. This branch remains a draft until the mathematical and 3D validation gates pass. **Publication readiness: 5/10, unchanged.**
