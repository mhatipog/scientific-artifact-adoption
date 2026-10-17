"""Offline checks for repeated-seed summary, not original model accuracy."""
import importlib.util
from pathlib import Path
import pytest

source = Path(__file__).resolve().parents[1] / 'seed_analysis.py'
spec = importlib.util.spec_from_file_location('seed_analysis', source)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


def trial(seed, correct=65, n=100):
    return {'seed': seed, 'correct': correct, 'n': n, 'accuracy': correct / n}


def test_three_seed_baseline_values_and_threshold_status():
    trials = [trial(0, 689, 1049), trial(1, 689, 1049), trial(2, 677, 1049)]
    v = s.summarize_trials(trials, 0.65)
    assert v['mean_accuracy'] == pytest.approx(0.6530028598665396)
    assert v['samples_per_seed'] == 1049
    assert v['seeds_at_or_above_historical_floor'] == 2
    assert v['seeds_below_historical_floor'] == 1
    assert v['mean_vs_historical_floor'] == 'at_or_above'
    assert v['standard_deviation_across_seeds'] > 0
    assert v['mean_accuracy_95pct_monte_carlo_t_interval'][0] < v['mean_accuracy'] < v['mean_accuracy_95pct_monte_carlo_t_interval'][1]


def test_all_below_floor_produces_completed_analytical_summary():
    v = s.summarize_trials([trial(i, 60) for i in range(30)], 0.65)
    assert v['mean_vs_historical_floor'] == 'below'
    assert v['seeds_below_historical_floor'] == 30
    assert v['standard_deviation_across_seeds'] == 0
    assert v['mean_accuracy_95pct_monte_carlo_t_interval'] == pytest.approx([0.6, 0.6])


def test_one_seed_has_no_mean_interval():
    v = s.summarize_trials([trial(0, 70)], 0.65)
    assert v['mean_accuracy_95pct_monte_carlo_t_interval'] is None
    assert v['seeds_at_or_above_historical_floor'] == 1


def test_unscored_smoke_summary():
    v = s.summarize_trials([trial(0, 70), trial(1, 72)], None)
    assert v['historical_single_run_floor'] is None
    assert 'seeds_at_or_above_historical_floor' not in v


@pytest.mark.parametrize('trials', [
    [], [trial(0), trial(0)], [trial(0, 70, 100), trial(1, 70, 101)],
    [dict(trial(0), accuracy=float('nan'))],
    [dict(trial(0), accuracy=1.4)],
])
def test_malformed_trials_rejected(trials):
    with pytest.raises(ValueError):
        s.summarize_trials(trials, 0.65)


def test_mean_interval_gets_narrower_with_more_independent_seed_evaluations():
    small = [trial(i, c) for i, c in enumerate([60, 70, 62, 68])]
    large = [trial(i, c) for i, c in enumerate([60, 70, 62, 68] * 10)]
    ci_small = s.summarize_trials(small, 0.65)['mean_accuracy_95pct_monte_carlo_t_interval']
    ci_large = s.summarize_trials(large, 0.65)['mean_accuracy_95pct_monte_carlo_t_interval']
    assert ci_large[1] - ci_large[0] < ci_small[1] - ci_small[0]
