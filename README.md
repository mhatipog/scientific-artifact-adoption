# Scientific artifact adoption: BreizhCrops reproduction

An independent, **unofficial** scientific reproducibility project. Its first target is the [BreizhCrops crop-classification benchmark](https://github.com/dl4sits/BreizhCrops) associated with Rußwurm et al., *BreizhCrops: A Time Series Dataset for Crop Type Mapping* (ISPRS Archives, 2020), [DOI 10.5194/isprs-archives-XLIII-B2-2020-1545-2020](https://doi.org/10.5194/isprs-archives-XLIII-B2-2020-1545-2020).

## Current evidence

- The upstream Transformer source at commit `6de796ed36a457c8520322d6110b8f2862fd8c25` was tested locally with modern CPU PyTorch, with correctly shaped, normalized outputs.
- This **does not** establish published classification accuracy. Full-dataset, original-checkpoint reproduction remains to be verified.
- A manual GitHub Actions workflow performs the authentic Belle-Île test and stores small JSON evidence reports. It does not redistribute any original source dataset or model checkpoint.

## Run

Open **Actions → BreizhCrops scientific reproduction → Run workflow**. The `run_real_benchmark` input defaults to `false`. Choose `true` only to opt in to original S3 downloads and real-data evaluation on an ephemeral CPU runner.

The workflow:
1. Clones original BreizhCrops at the exact historical commit (not a rewritten model).
2. Installs pinned Python dependencies.
3. Downloads the *original author's* checkpoint using the original recorded URL.
4. Loads 2017 L1C Belle-Île using the upstream dataset loader.
5. Runs three seeded, full-dataset evaluations and writes checkpoint/dataset hashes and accuracies to `workflow-results/breizhcrops-real.json`.
6. Uploads JSON evidence only, whether the attempt passes or encounters a blocker.

The upstream test in `tests/models.py` expects a Transformer Belle-Île accuracy >= **0.69 - 0.04 = 0.65**. The report distinguishes `pass`, `below_historical_threshold`, and `blocked`. The original transform randomly subsamples temporal observations, so reported accuracy may vary with the seed.

### Interpretation

`pass` means this particular original-weight, original-data evaluation reached the historical regression threshold on all recorded seeds. It is **not** a claim to have reproduced every experiment in the paper. `blocked` and synthetic-only tests never count as experimental replication.

The benchmark uses GitHub-hosted internet access and consumes external download bandwidth. The upstream authors control the checkpoint and dataset URLs; this project does not mirror those files. No automated scheduled jobs or GPU runs.

## Attribution and licensing

Original data, code, checkpoints, paper and their associated license remain with their authors. This independent harness does not claim authorship of BreizhCrops and does not relicense their work. For exact model and preprocessing behavior, follow the historical source revision linked above.

Related follow-on research on thin-strip contact-line solvers and bendotaxis is kept separately until there is independently verified upstream replication evidence. The project aims to submit focused upstream fixes after the evidence is collected.

## Upstream integration

Possible downstream adoption: [TorchGeo request #3618](https://github.com/torchgeo/torchgeo/issues/3618) (BreizhCrops dataset). Any future integration should preserve the original benchmark's dataset semantics and attribution.
