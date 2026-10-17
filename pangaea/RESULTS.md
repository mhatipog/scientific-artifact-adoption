# PANGAEA evaluator invariance audit — 17 October 2026

**Reproduced against the real upstream source**, pinned to Git commit
`ca19b837a16797cd809225225b671473c7c93768` and evaluator blob
`345e7b33a471134e2d1a5bd3f1d28c9d5776b73c`.

[Successful GitHub Actions experiment and downloadable evidence](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39035890791).

The workflow clones the original PANGAEA source, checks both Git hashes,
extracts the original evaluator methods *unchanged* via AST and runs their
metric calculations using real PyTorch operations. External data loading
and heavy model dependencies are replaced by fixed synthetic predictions.
It does **not** recreate any trained foundation-model checkpoint.

## Findings

**1. Regression:** `RegEvaluator.evaluate` averages batch MSEs separately,
all-reduces their sum across ranks, and divides by the local batch count.
Consequently nine fixed predictions from A with true MSE 1.00 report
MSE 4.50 with batch size 8, while fixed predictions from B with true MSE
2.00 report 2.00. Ranking reverses *solely* because of batching.
When eight fixed squared errors are split across two real CPU/Gloo ranks,
the original code reports MSE **9.00** rather than the true **4.50**.
A 39-line guarded correction accumulates global squared-error **sum and
pixel count** instead; it passes `git apply --check` against the pinned
original source and reproduces the mathematically correct scores.

**2. Segmentation:** Original confusion matrix has prediction rows and
reference columns. `SegEvaluator.compute_metrics` consequently labels
per-class recall as **precision** and precision as **recall**. For
`[[40,10],[20,30]]`, upstream class-0 `Precision` is 66.67%
(the actual recall) whereas the true precision is **80%**. The matrix
has been checked directly with scikit-learn. The class-label mix-up does
not itself alter mIoU or F1.

**3. Multi-label KNN:** `KNNClassificationEvaluator.evaluate` computes
micro-F1 separately on each batch and averages these batch scores.
Four fixed multi-label predictions produce reported micro-F1 **0.5000**
with batch size 1 and **0.6667** with batch size 4. The actual global
micro-F1 is **0.6667**, independently verified with scikit-learn.

## Scientific impact and essential cautions

PANGAEA is a contemporary benchmark for geospatial foundation models,
including forest-biomass regression. The affected evaluator Git blob
also appears in the NASA ILAB and AGBD downstream forks, confirming
**source inheritance**, not actual execution of the affected path.

- PANGAEA defaults to `test_batch_size: 1`. Under single-rank evaluation
  this avoids the underfilled-batch issue. The multi-rank multiplication
  remains possible, but a fixed scaling factor alone does **not**
  necessarily reverse rankings if every model uses identical execution.
- The 2026 AGBD fork's experiment scripts use **`AGBDEvaluator`**, a
  separate class that correctly globally sums errors and counts. Its
  published results must **not** be alleged erroneous based on merely
  sharing the older base evaluator file.
- We have not rerun a trained PANGAEA checkpoint or official BioMassters
  predictions. **No published paper or leaderboard is claimed incorrect.**
- Padded/duplicated distributed dataset sampling and invalid-value masks
  require separate auditing; this minimal patch does not solve them.

The contribution is a **benchmark metric-invariance contract**: frozen
predictions should produce the same global metric under equivalent batch
partition or distributed execution, with unambiguous class-axis semantics.

Next action: offer the regression fix to the original maintainers in a
concise issue. Treat the other two corrections separately. Then examine
genuinely released model predictions and splits before any publication claim.
