# Run the first authentic BreizhCrops experiment

Repository: https://github.com/mhatipog/scientific-artifact-adoption

1. Open https://github.com/mhatipog/scientific-artifact-adoption/actions/workflows/breizhcrops-reproduction.yml
2. Select **Run workflow**, initially leave **run_real_benchmark = false**. The first run tests original-model CPU compatibility on synthetic inputs only.
3. Download the **breizhcrops-reproduction-evidence** artifact and confirm the synthetic report.
4. Run again with **run_real_benchmark = true**. This explicitly opts in to downloading the published checkpoint and original 2017 Belle-Île data.
5. Download the same artifact, read **breizhcrops-real.json**. Successful historical regression evidence requires status=pass, three full-dataset seeds, max_samples=null, original source revision, asset SHA-256 hashes, and accuracy >=0.65 for each seed.

If blocked, retain the failure log and resolve that precise issue. Never replace original assets to force a green result. No scheduling, GPU, or dataset uploads. A 25-minute CPU cap applies. Actions usage depends on your account settings.

When we have real evidence, propose a focused upstream change at https://github.com/dl4sits/BreizhCrops, then coordinate TorchGeo https://github.com/torchgeo/torchgeo/issues/3618.
