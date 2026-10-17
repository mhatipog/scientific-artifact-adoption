"""Offline tests of the harness; NOT published-model/data reproduction."""
import importlib.util
from pathlib import Path
import json
import numpy as np
import torch

SOURCE = Path(__file__).resolve().parents[1] / "run_benchmark.py"
spec = importlib.util.spec_from_file_location("runner", SOURCE)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)

class FakeDataset(torch.utils.data.Dataset):
    def __len__(self):
        return 7
    def __getitem__(self, idx):
        values = np.random.choice(np.arange(10), 45, replace=True)
        x = torch.tensor(np.broadcast_to(values[:, None], (45, 13)).copy(), dtype=torch.float32)
        return x, torch.tensor(idx % 9, dtype=torch.long), idx

class FakeModel(torch.nn.Module):
    def forward(self, x):
        scores = torch.zeros((len(x), 9))
        scores[:, 0] = x[:, 0, 0]
        return scores.log_softmax(dim=-1)

def test_repeatability_and_shapes():
    a = r.evaluate(FakeModel(), FakeDataset(), seed=11, batch_size=3, max_samples=None)
    b = r.evaluate(FakeModel(), FakeDataset(), seed=11, batch_size=3, max_samples=None)
    assert a == b and a["n"] == 7 and 0 <= a["accuracy"] <= 1

def test_subset_is_not_full_benchmark():
    a = r.evaluate(FakeModel(), FakeDataset(), seed=0, batch_size=2, max_samples=2)
    assert a["n"] == 2

def test_missing_checkpoint_is_blocked(tmp_path):
    import argparse
    args = argparse.Namespace(upstream_checkout=None, dataset_root=tmp_path/"data",
        region="belle-ile", weights=tmp_path/"missing.pth",
        download_official_weights=False, offline=True,
        batch_size=2, seeds=[0], max_samples=None, output=tmp_path/"report.json",
        require_upstream_revision=False)
    assert r.run(args) == 1
    report = json.loads(args.output.read_text())
    assert report["status"] == "blocked" and report["trials"] == []
    assert "mean_accuracy" not in report

def test_historical_thresholds_are_from_original_tests():
    assert r.HISTORICAL_MINIMUM[("belle-ile", "L1C")] == 0.65
    assert r.HISTORICAL_MINIMUM[("frh04", "L1C")] == 0.78
