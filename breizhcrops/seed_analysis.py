"""Analysis of repeat inference stochasticity on one fixed pretrained model/dataset.

This is NOT a confidence interval for population generalization or model training.
Each seed drives the original upstream random temporal subsampling for the
same fixed parcels, labels, checkpoint, and model.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import t as student_t


def summarize_trials(trials: list[dict], historical_minimum: float | None) -> dict:
    """Produce reproducible descriptive statistics and a per-seed t interval.

    The t interval is an approximate Monte Carlo interval around mean accuracy
    under the source dataset's random temporal sampling procedure.
    """
    if not trials:
        raise ValueError("At least one completed trial is required")
    seeds = [int(t["seed"]) for t in trials]
    if len(set(seeds)) != len(seeds):
        raise ValueError("Repeated random seed")
    sizes = [int(t["n"]) for t in trials]
    if len(set(sizes)) != 1 or sizes[0] <= 0:
        raise ValueError("Trials must cover the same positive number of samples")
    accuracies = np.asarray([float(t["accuracy"]) for t in trials], dtype=float)
    if not np.isfinite(accuracies).all() or np.any((accuracies < 0) | (accuracies > 1)):
        raise ValueError("Accuracy must be a finite probability")
    mean = float(accuracies.mean())
    standard_deviation = (0.0 if np.all(accuracies == accuracies[0]) else float(accuracies.std(ddof=1))) if len(trials) > 1 else None
    ci = None
    if standard_deviation is not None:
        half_width = float(student_t.ppf(0.975, len(trials) - 1) * standard_deviation / np.sqrt(len(trials)))
        ci = [max(0.0, mean - half_width), min(1.0, mean + half_width)]
    q05, q25, q50, q75, q95 = np.quantile(accuracies, [0.05, 0.25, 0.50, 0.75, 0.95]).tolist()
    result = {
        "number_of_seeds": len(trials),
        "samples_per_seed": sizes[0],
        "mean_accuracy": mean,
        "min_accuracy": float(accuracies.min()),
        "max_accuracy": float(accuracies.max()),
        "standard_deviation_across_seeds": standard_deviation,
        "mean_accuracy_95pct_monte_carlo_t_interval": ci,
        "accuracy_quantiles": {"p05": q05, "p25": q25, "p50": q50, "p75": q75, "p95": q95},
        "historical_single_run_floor": historical_minimum,
        "interpretation": (
            "Fixed dataset, labels, and official checkpoint; randomness comes from "
            "upstream temporal subsampling only. The across-seed interval estimates "
            "Monte Carlo uncertainty in mean accuracy, NOT uncertainty from field "
            "sampling, model fitting, or deployment to another geography."
        ),
    }
    if historical_minimum is not None:
        if not (0 <= historical_minimum <= 1):
            raise ValueError("Historical threshold outside [0, 1]")
        passed = int(np.sum(accuracies >= historical_minimum))
        result.update({
            "seeds_at_or_above_historical_floor": passed,
            "seeds_below_historical_floor": len(trials) - passed,
            "fraction_seeds_at_or_above_historical_floor": passed / len(trials),
            "mean_vs_historical_floor": "at_or_above" if mean >= historical_minimum else "below",
            "historical_floor_is_stochastic_single_run_guardrail": True,
        })
    return result
