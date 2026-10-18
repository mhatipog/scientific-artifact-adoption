# Guarded fifth-order phase design (Round 14)

This is a **theoretical inverse-design audit** in the existing contact-angle hysteresis programme. It is **not** part of the experimental Yan image-volume reconstruction. The latter remains separately documented in draft research PR #14.

- [Full scientific note: mathematical, numerical and novelty qualifications](../ROUND14_GUARDED_PHASE_PAIR.md)
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

## Exact single-front pressure–volume bridge (Round 15)

The [Round-15 theorem, proof, numerical cross-check and conditional pressure-noise certificate](../ROUND15_SINGLE_FRONT_BRIDGE.md) show that at a stationary pressure, within this reduced two-front model,

\[
B_4=-[H'(s_c)]^4\,P'''(Q_c)/12,\qquad Q=H(s),\;P(Q)=p(s).
\]

This eliminates explicit \(H''\) dependence **exactly**. The sign of the single-front *third pressure–volume derivative* fixes the reduced local two-front quartic class if the other modes are stable. A [pure-standard-library implementation](pressure_volume_bridge.py) checks the chain rule and provides a conservative four-pressure-point sign bound, with [five zero-network regression tests](tests/test_pressure_volume_bridge.py). No independent calibrated \(P(Q)\) data, nonneutral 3D proof or physical snap prediction is claimed.

Run offline:

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_pressure_volume_bridge.py' -v
```
