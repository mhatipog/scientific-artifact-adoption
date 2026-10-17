# Proposed upstream report: fast Transformer test is stochastic

**Title:** Fast Transformer test varies with random temporal sampling (8/30 seeds below threshold)

The default BreizhCrops evaluation transform samples 45 observations with `np.random.choice`, and `tests/models.py` checks Transformer Belle-Île accuracy against **0.65** without fixing a seed.

I ran the original 2022 source and official checkpoint on **1,049** Belle-Île fields (2017 L1C) with PyTorch 2.10 CPU, changing only the global random seed (0–29). Results: **22/30 passes**, **8/30 below threshold**, accuracy **62.73–66.54%** (mean **65.14%**). Original model, dataset loader, and temporal transform were unchanged.

Evidence: [all 30 seed-by-seed counts and asset hashes](https://github.com/mhatipog/scientific-artifact-adoption/blob/main/evidence/2026-10-17-belle-ile-30-seeds.json); [successful CPU workflow](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39021581078).

This suggests a potentially flaky test because its temporal subsampling is stochastic. I haven't compared against the original Python/PyTorch environment, so absolute-accuracy differences may have another cause.

Would you prefer fixed deterministic sampling for the fast test, or a multi-seed criterion? Happy to submit a small targeted PR.

---

**Status:** Draft only. The GitHub app has no permission to create an issue in the original research repository. Nothing has been posted upstream.
