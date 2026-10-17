# Scientific artifact adoption: BreizhCrops reproduction

An independent, **unofficial** scientific reproducibility project. Its first target is the [BreizhCrops crop-classification benchmark](https://github.com/dl4sits/BreizhCrops) associated with Rußwurm et al., *BreizhCrops: A Time Series Dataset for Crop Type Mapping* (ISPRS Archives, 2020), [DOI 10.5194/isprs-archives-XLIII-B2-2020-1545-2020](https://doi.org/10.5194/isprs-archives-XLIII-B2-2020-1545-2020).

## Repeated-seed result (17 October 2026)

**Completed real-data reproduction experiment**: the original pretrained Transformer, upstream source revision `6de796e`, and 2017 L1C Belle-Île dataset (1,049 fields) ran successfully in a CPU GitHub Actions job over 30 fixed seeds. No original model/dataset/transform changes.

| Metric | Observed |
| --- | ---: |
| Mean accuracy across random temporal samplings | **65.14%** |
| Standard deviation across seeds | 0.82 percentage points |
| Minimum / maximum seed accuracy | 62.73% / 66.54% |
| Seeds at or above the historical 65% single-run floor | **22 / 30** |
| Seeds below the floor | **8 / 30** |
| Approximate 95% Monte Carlo interval for mean | 64.84%–65.45% |

The 65% floor comes from the upstream fast Transformer test's `0.69 - 0.04` regression tolerance. It is **not a guarantee across all seeds**. Its stochastic pass/fail outcome is meaningful for regression-test design. The interval only covers random temporal subsampling on the same fields and checkpoint, **not generalization to other field datasets**. The reported mean does not demonstrate the original nominal 69% expected value.

[Permanent 30-seed results](evidence/2026-10-17-belle-ile-30-seeds.json) · [successful workflow run](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39021581078) · [original 3-seed run](https://github.com/mhatipog/scientific-artifact-adoption/actions/runs/39020692334)

**The workflow has been returned to manual-only mode.** No scheduled or push-triggered jobs are left active.

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

## Re-run the experiment

For an independent rerun, navigate to [GitHub Actions](https://github.com/mhatipog/scientific-artifact-adoption/actions/workflows/breizhcrops-reproduction.yml), click **Run workflow**, set `run_real_benchmark` to **true**, and select `main`. The full test uses 30 deterministic seeds of the original stochastic evaluation procedure, and uploads JSON evidence only. No scheduled or push triggers remain. The evidence artifact expires after 14 days; the summarized original run is preserved permanently under `evidence/`.


## Contact-angle hysteresis: two-front statics (2026)

The independent [CAH reduced-order reproducibility study](cah-statics/README.md)
extends the 2019 microchannel statics theory to two coupled moving menisci.
It reproduces an internally checked example and publication-quality
figures, but **does not yet establish quantitative agreement with the
published experiments**. Run the self-contained manual/pr-scoped
[CAH statics Actions workflow](.github/workflows/cah-statics.yml) for the
outputs. This research is separate from the BreizhCrops artifact adoption.
