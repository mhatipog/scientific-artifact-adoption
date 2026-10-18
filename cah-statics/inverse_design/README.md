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

## Exact single-front pressure–volume bridge (Round 15)

The [Round-15 independent rederivation, numerical cross-check and conditional pressure-noise certificate](../ROUND15_SINGLE_FRONT_BRIDGE.md) show that at a stationary pressure, within this reduced two-front model,

\[
B_4=-[H'(s_c)]^4\,P'''(Q_c)/12,\qquad Q=H(s),\;P(Q)=p(s).
\]

The identity was already documented in [the earlier unmerged pressure-volume transfer note](https://github.com/mhatipog/scientific-artifact-adoption/blob/research/cah-phase-volume-transfer-20261018/cah-statics/theory/PRESSURE_VOLUME_TRANSFER.md); this round contributes a separate implementation, the more resolvable geometry cross-check and the conditional measurement-error bound. This eliminates explicit \(H''\) dependence **exactly**. The sign of the single-front *third pressure–volume derivative* fixes the reduced local two-front quartic class if the other modes are stable. A [pure-standard-library implementation](pressure_volume_bridge.py) checks the chain rule and provides a conservative four-pressure-point sign bound, with [five zero-network regression tests](tests/test_pressure_volume_bridge.py). No independent calibrated \(P(Q)\) data, nonneutral 3D proof or physical snap prediction is claimed.

Run offline:

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_pressure_volume_bridge.py' -v
```

## Direct nonlinear branch check

The [Round-15 note](../ROUND15_SINGLE_FRONT_BRIDGE.md) also compares the measured \(\eta^2/(q-Q_c)\) of eight **directly solved untruncated reduced-model equilibria** against \(-6P''/P'''\). The relative error falls to **0.065% / 0.022%** for a 0.5-µm centre shift in the −/+ geometries, respectively, with \(O(c^2)\) convergence. The [solver and quadrature reproducer](local_branch_validation.py), [eight-row numerical data](derived/nonlinear_local_branch.csv) and [three offline regression tests](tests/test_local_branch_validation.py) contain the exact outputs. This tests the *reduced local equilibrium branch*, not its 3D stability or a jump distance.

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_local_branch_validation.py' -v
python cah-statics/inverse_design/local_branch_validation.py --csv
```

## Neutral-plate continuum spectral isolation (Round 16)

The [Round-16 mathematical audit](ROUND16_SPECTRAL_ISOLATION.md) establishes a
**conditional finite-spectral non-identifiability result**: signed Fourier
roughness phases with exactly matched finite-order correlations can produce
opposite *nonlinear* local quartics even as every *fixed* finite collection
of free-contact-line **neutral-plate 3D Jacobi eigenvalues** approaches
agreement. It explicitly includes the solid wall's support curvature in
the Robin contact-line condition; using the critical-point Robin coefficient
away from a pressure stationary state would give spurious zero modes.

The nominal 170-nm 5/7 pair's first positive eigenvalues are
12.526282 / 12.509126 mm^-2, and the first depth-varying mode at
d=0.38 mm is 68.349061 mm^-2. The independent P1 FEM gaps show
second-order convergence. These are **exact special-case continuum linear
spectra, not finite-angle nonlinear 3D quartics**. The theorem is a
synthesis of our prior Round-8 Fourier argument with conventional
Sturm–Liouville theory; literature priority has not been established.

- [Executable continuum spectral calculation and CSV reproducer](spectral_isolation.py)
- [Ten offline FEM, boundary-condition and asymptotic tests](tests/test_spectral_isolation.py)
- [Derived nominal and high-harmonic CSV](derived/neutral_plate_jacobi_spectrum.csv)
- [Reproducible SVG figure](figures/jacobi_spectral_isolation.svg) and [generator](plot_spectral_isolation.py)

Run locally, without JAX, GitHub Actions or experimental-media downloads:

```sh
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_spectral_isolation.py' -v
python cah-statics/inverse_design/spectral_isolation.py --csv
python cah-statics/inverse_design/plot_spectral_isolation.py
```

**Still draft and unmerged.** No fabricated pair, measured wetting or
non-neutral constrained nonlinear 3D stability. Submission readiness remains 5/10.

## Exact neutral-plate plate-wetting response and invisible quartic sign (Round 17)

[Round-17 analytical note, convergent FEM, constructive jet counterexample, and failed nonlinear gate](ROUND17_NEUTRAL_PLATE_RESPONSE.md):

- For the selected model cylinders, the fixed-volume plate-angle energy coefficient is **Q− = −0.087315948234** and **Q+ = −0.085192444372 mm²**, with rigorous spectral truncation bounds below **1.85×10⁻¹¹ mm²** for 320 modes (aside from floating and input errors). Independent bordered P1 FEM shows approximately fourfold error reduction per mesh refinement; 11 new offline tests pass.
- A mirror-even smooth fourth-jet wall perturbation changes **B4** by exactly **H′sin(α) δw⁽⁴⁾/[12w(1+(w′/2)²)]** without changing the selected neutral-plate Jacobi operator, its entire spectrum or **Q**. This is a mathematical inverse-identifiability counterexample, **not** a realizable fabrication design or a new generic normal-form principle; seven more offline tests pass.
- A separately implemented full-area CasADi nonlinear 3D model **fails to converge to the exact 90° corrugated pressure**, despite very small KKT residuals. Do not treat its 89° energy coefficients or bifurcation classes as validated. [Rejected mesh diagnostic](derived/nonlinear3d_failed_mesh_gate.csv), [optional reproducer](nonlinear3d_casadi_audit.py).

Reproduce mathematics without workflows or JAX:

```sh
python cah-statics/inverse_design/plate_susceptibility.py --csv
python cah-statics/inverse_design/jet_invisibility.py
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_plate_susceptibility.py' -v
python -m unittest discover -s cah-statics/inverse_design/tests -p 'test_jet_invisibility.py' -v
```

**Still a draft:** no nonneutral continuum nonlinear quartic, independent wettability measurement, fabricated pair or optical liquid-depth calibration; readiness **5/10, unchanged**.

## Geometry-fitted nonlinear 3D continuation (Round 18)

[Full audited note](ROUND18_GEOMETRY_FITTED_3D.md), [reproducible solver](nonlinear3d_local_continuation.py), [seven offline tests](tests/test_nonlinear3d_local_continuation.py), [raw 10-grid convergence CSV](derived/nonlinear3d_local_diagonal.csv), [independent axial/depth grid CSV](derived/nonlinear3d_local_grid_axes.csv), [two-axis Richardson CSV](derived/nonlinear3d_local_richardson.csv) and [generated SVG figure](figures/nonlinear3d_local_convergence.svg).

Replacing uniform transverse sampling by nodes graded in the **physical axial coordinate of the exact neutral cylinder** and using a **local sparse KKT equilibrium solve** resolves the earlier numerical branch/aliasing difficulty for near-neutral 3D single menisci. At 89.9° plates, finite-grid nonlinear Q differs by **0.56%/0.45%** from the independent exact continuum Jacobi susceptibilities; exploratory double-direction Richardson errors are **0.016%/0.006%**. Exact 90° pressure residuals are ≲1×10⁻⁵ mm⁻¹ on the finest 257×17 mesh. Discrete constrained Hessian signs are positive on separately tested 65×5 and 129×9 neutral/nonneutral states.

**These are numerical benchmark results, not a measured effect, continuum stability theorem or nonneutral two-front B4.** The existing failed unconstrained 3D result remains in Round 17 as an essential negative control. No workflows, media or experimental pressure fits were introduced. **Still draft/unmerged; journal readiness 5/10, unchanged.**
