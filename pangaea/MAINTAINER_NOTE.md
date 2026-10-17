# Draft for upstream PANGAEA issue — not yet posted

**Title:** Regression MSE depends on batch size and distributed world size

I was reproducing the regression evaluator on fixed predictions and noticed
that it averages each batch separately, then sums across ranks but divides
by the local number of batches.

For nine examples the original evaluator reports MSE 4.5 rather than 1.0
when the final batch contains one example. With two CPU/Gloo ranks,
it reports 9.0 instead of a global MSE of 4.5.

I checked this against `RegEvaluator.evaluate` in the original pinned
source and verified that accumulating squared errors and pixel counts
across ranks resolves both:

[Reproduction, 39-line patch and successful CI](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39035890791)

Would a small PR for the regression accumulator be useful?

No inference about the paper's reported results is intended.
