# Three-dimensional Young–Laplace model — reproducible benchmark

See [Round 5 research memorandum](ROUND5_3D_NONNEUTRAL.md) for detailed derivation, assumptions, and publishability review.

This code is a novel finite-depth Galerkin model based on four-wall capillary surface energy. Results are **model-only**, not experimentally validated.

Dependencies: `numpy scipy matplotlib jax[cpu]` with Python >=3.11.

Use `python surface_3d.py`, `python full_stability_3d.py` and `python figures_round5.py`. The continuation dataset CSV is included. The 60–70 degree threshold is **not converged and not a validated finding**.