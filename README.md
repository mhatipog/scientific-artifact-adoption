# Scientific artifact adoption: BreizhCrops reproduction

An independent, **unofficial** scientific reproducibility project. Its first target is the [BreizhCrops crop-classification benchmark](https://github.com/dl4sits/BreizhCrops) associated with Rußwurm et al., *BreizhCrops: A Time Series Dataset for Crop Type Mapping* (ISPRS Archives, 2020), [DOI 10.5194/isprs-archives-XLIII-B2-2020-1545-2020](https://doi.org/10.5194/isprs-archives-XLIII-B2-2020-1545-2020).

## Current evidence

- The upstream Transformer source at commit `6de796ed36a457c8520322d6110b8f2862fd8c25` was tested locally with modern CPU PyTorch, with correctly shaped, normalized outputs.
- On 17 October 2026, an actual GitHub Actions run evaluated the original pretrained checkpoint on all **1,049** Belle-Île fields, using the original data and source code. Across seeds 0/1/2 it measured **65.68%, 65.68%, 64.54%** (mean **65.30%**). The first workflow returned failure because it incorrectly required *every seed* to meet the upstream **single-run** 65% minimum. [Preserved baseline evidence](evidence/2026-10-17-belle-ile-three-seed-baseline.json) · [Original Actions run](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39020692334).
- A manual GitHub Actions workflow performs the authentic Belle-Île test and stores small JSON evidence reports. It does not redistribute any original source dataset or model checkpoint.

## Run

Open **Actions → BreizhCrops scientific reproduction → Run workflow**. The `run_real_benchmark` input defaults to `false`. Choose `true` only to opt in to original S3 downloads and real-data evaluation on an ephemeral CPU runner.

The workflow:
1. Clones original BreizhCrops at the exact historical commit (not a rewritten model).
2. Installs pinned Python dependencies.
3. Downloads the *original author's* checkpoint using the original recorded URL.
4. Loads 2017 L1C Belle-Île using the upstream dataset loader.
5. Runs **30 fixed-seed**, full-dataset evaluations (seeds 0–29) and records all per-seed accuracies, dataset/checkpoint hashes, descriptive variation, and an approximate 95% Monte Carlo t interval around the mean in `workflow-results/breizhcrops-real.json`.
6. Uploads JSON evidence only, whether the attempt passes or encounters a blocker.

The upstream test in `tests/models.py` checks a **single random inference run** against **0.69 - 0.04 = 0.65**. The original transform resamples temporal observations, so one seed can fall below 65% while others exceed it. The new report distinguishes `completed` (experiment ran, regardless of accuracy), `smoke_only`, and `blocked`. It separately records **mean accuracy relative to the historical floor, standard deviation, empirical quantiles, and the proportion of seeds clearing the floor**. Valid below-floor results must never be silently converted into passing scientific results.

### Interpretation

`completed` means authentic data and the official checkpoint were evaluated for every requested seed. A successful workflow means the experiment **executed correctly**, not that every seed (or even the mean) reached the historical accuracy. The scientific comparison is recorded in `seed_analysis.mean_vs_historical_floor` and per-seed results. A 95% t interval summarizes Monte Carlo uncertainty from randomized temporal subsampling **on the same 1,049 fields**; it is **not** an interval for generalization accuracy on new farms, countries, training sets, or distributions. `blocked` and synthetic-only tests never count as real-data replication.

The benchmark uses GitHub-hosted internet access and consumes external download bandwidth. The upstream authors control the checkpoint and dataset URLs; this project does not mirror those files. No automated scheduled jobs or GPU runs.

## Attribution and licensing

Original data, code, checkpoints, paper and their associated license remain with their authors. This independent harness does not claim authorship of BreizhCrops and does not relicense their work. For exact model and preprocessing behavior, follow the historical source revision linked above.

Related follow-on research on thin-strip contact-line solvers and bendotaxis is kept separately until there is independently verified upstream replication evidence. The project aims to submit focused upstream fixes after the evidence is collected.

## Upstream integration

Possible downstream adoption: [TorchGeo request #3618](https://github.com/torchgeo/torchgeo/issues/3618) (BreizhCrops dataset). Any future integration should preserve the original benchmark's dataset semantics and attribution.

## Next experiment

The manually dispatched real-data workflow now runs 30 fixed seeds. Navigate to [Actions](https://github.com/mhatipog/scientific-artifact-adoption/actions/workflows/breizhcrops-reproduction.yml), click **Run workflow**, set `run_real_benchmark` to **true** and use branch `main`. No schedule or push trigger is enabled, so expensive downloads cannot start unexpectedly. Retrieve the JSON evidence artifact; it is retained for 14 days. Do not interpret 30 Monte Carlo seeds as 30 independent geographic test sets.
