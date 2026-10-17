# Reproduction of PANGAEA regression evaluation invariance

This is a **source-level verification**, deliberately independent from model
training and remote-sensing data downloads.

PANGAEA's original *RegEvaluator.evaluate* is extracted unchanged by AST
from the **actual pinned upstream source** (not a reimplementation).
The original and a narrowly modified version execute on the same fixed
predictions, with actual PyTorch/Gloo collective operations.

Original source: https://github.com/VMarsocci/pangaea-bench
Revision: `ca19b837a16797cd809225225b671473c7c93768`
Blob: `345e7b33a471134e2d1a5bd3f1d28c9d5776b73c`

## Invariants under test

1. Evaluation batch size cannot change MSE for frozen predictions.
2. Evaluation world size cannot change MSE for the same complete collection
   of frozen predictions.
3. A model comparison should be unaffected by the grouping of samples.

The simple reproducible counterexample has 9 scalar predictions:
- A squared errors: 8 zeros and one 9, global MSE=1.
- B squared errors: nine 2s, global MSE=2.

With batch size 8, upstream computes mean([0,9]) = 4.5 for A,
so it reverses the true A-before-B ranking. With batch size 1 the
upstream calculation gives the correct global MSE in one process.

With 2 ranks, 4 equally sized samples/rank, squared errors 1..8,
upstream returns 9.0; correct global MSE is 4.5.

Only the three accumulator code spans change in the corrected version:
store squared-error **sum** and pixel **count** as FP64 totals,
all-reduce both, then divide their global totals. No image data, labels
or model are changed. The correction also guards the empty case.

This project has NO benchmark checkpoint or frozen real predictions;
it cannot establish that published PANGAEA model rankings were
incorrect. Note that if all model candidates use the same rank count
and equal batches, the world-size multiplier alone may not alter their
order, even though the reported absolute MSE/RMSE are wrong.

Run: inspect the workflow, then execute it via push or workflow_dispatch
on this isolated research branch in scientific-artifact-adoption.

Next: submit a minimal maintainer-facing proposal once verification
is complete; separately pursue a real checkpoint + data replication.
