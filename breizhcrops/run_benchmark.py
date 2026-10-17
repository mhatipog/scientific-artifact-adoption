#!/usr/bin/env python3
"""Reproducibility runner for the historical BreizhCrops Transformer evaluation.

Independent harness: imports original BreizhCrops code from a pinned checkout.
Does not modify original source, model, preprocessing, or training protocol.
Scientific provenance: Russwurm et al. ISPRS Archives 2020, DOI 10.5194/isprs-archives-XLIII-B2-2020-1545-2020.
Maintainer-facing tool, not a claim of completed replication.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import urllib.request
import numpy as np
import torch
from torch.utils.data import DataLoader

UPSTREAM_SHA = "6de796ed36a457c8520322d6110b8f2862fd8c25"
DOI = "10.5194/isprs-archives-XLIII-B2-2020-1545-2020"
OFFICIAL_WEIGHTS_URL = (
    "https://breizhcrops.s3.eu-central-1.amazonaws.com/models/eval/TransformerEncoder/"
    "TransformerEncoder_input-dim%3D13_num-classes%3D9_d-model%3D64_d-inner%3D128_"
    "n-layers%3D5_n-head%3D2_dropout%3D0.017998950510888446_"
    "learning-rate%3D0.00017369201853408445_weight-decay%3D3.5156458637523697e-06.pth"
)
# Derived from upstream tests/models.py (expected 0.69 +/- 0.04 for Belle-Île).
HISTORICAL_MINIMUM = {("belle-ile", "L1C"): 0.65, ("frh04", "L1C"): 0.78}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temp, path)


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if hasattr(torch, "use_deterministic_algorithms"):
        torch.use_deterministic_algorithms(True, warn_only=True)


def git_revision(repo_path: Path | None):
    if repo_path is None:
        return None
    try:
        rev = subprocess.run(["git", "-C", str(repo_path), "rev-parse", "HEAD"], check=True,
                             capture_output=True, text=True, timeout=10).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(repo_path), "status", "--porcelain"], check=True,
                                    capture_output=True, text=True, timeout=10).stdout.strip())
        return {"revision": rev, "dirty": dirty, "matches_target": rev == UPSTREAM_SHA}
    except (subprocess.SubprocessError, OSError):
        return {"revision": "unknown", "dirty": None, "matches_target": None}


def download_official_weights(destination: Path):
    """Explicit opt-in only; no silent download or pickle unpickling."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return destination
    if not OFFICIAL_WEIGHTS_URL.startswith("https://breizhcrops.s3.eu-central-1.amazonaws.com/"):
        raise ValueError("Untrusted download endpoint")
    req = urllib.request.Request(OFFICIAL_WEIGHTS_URL, headers={"User-Agent": "breizhcrops-reproduction/0.1"})
    temp = destination.with_suffix(destination.suffix + ".part")
    try:
        with urllib.request.urlopen(req, timeout=45) as response, temp.open("wb") as output:
            if response.geturl().split("/")[2] != "breizhcrops.s3.eu-central-1.amazonaws.com":
                raise ValueError("Unexpected download redirect")
            total = 0
            while True:
                data = response.read(1024 * 1024)
                if not data:
                    break
                total += len(data)
                if total > 100 * 1024 * 1024:
                    raise ValueError("Unexpectedly large checkpoint (>100MB)")
                output.write(data)
        os.replace(temp, destination)
    finally:
        temp.unlink(missing_ok=True)
    return destination


def load_original_model(weight_path: Path):
    # Construct the actual upstream class, not a reimplementation.
    mod = importlib.import_module("breizhcrops.models")
    model = mod.TransformerModel(input_dim=13, num_classes=9).cpu()
    try:
        data = torch.load(weight_path, map_location="cpu", weights_only=True)
    except Exception as exc:
        raise RuntimeError("Checkpoint failed safe weights_only loading; no unsafe pickle fallback performed") from exc
    if not isinstance(data, dict):
        raise TypeError("Expected a dict checkpoint")
    weights = data.get("model_state", data)
    model.load_state_dict(weights, strict=True)
    model.eval()
    return model


def evaluate(model, dataset, *, seed: int, batch_size: int, max_samples: int | None):
    set_seed(seed)
    if max_samples is not None:
        if max_samples < 1:
            raise ValueError("max_samples must be positive")
        dataset = torch.utils.data.Subset(dataset, range(min(max_samples, len(dataset))))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)
    correct = 0
    seen = 0
    prediction_hasher = hashlib.sha256()
    labels = set()
    with torch.inference_mode():
        for x, y, field_id in loader:
            if not isinstance(x, torch.Tensor) or x.ndim != 3 or x.shape[1:] != (45, 13):
                raise ValueError(f"Expected [batch,45,13] after upstream transform, got {getattr(x,'shape',None)}")
            if not isinstance(y, torch.Tensor) or y.ndim != 1:
                raise ValueError("Unexpected labels from upstream dataset")
            logp = model(x.cpu())
            if not torch.isfinite(logp).all() or logp.shape != (len(y), 9):
                raise ValueError("Nonfinite or malformed model predictions")
            if not torch.allclose(logp.exp().sum(dim=1), torch.ones(len(y)), atol=1e-5):
                raise ValueError("Expected normalized log-probabilities")
            predicted = logp.argmax(dim=-1).cpu()
            truth = y.long().cpu()
            correct += int((predicted == truth).sum())
            seen += len(truth)
            labels.update(truth.tolist())
            for yi, pi in zip(truth.tolist(), predicted.tolist()):
                prediction_hasher.update(f"{yi}:{pi}\n".encode("ascii"))
    if not seen:
        raise ValueError("No samples evaluated")
    return {"seed": seed, "n": seen, "correct": correct, "accuracy": correct / seen,
            "labels_seen": sorted(labels), "prediction_sha256": prediction_hasher.hexdigest()}


def build_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--upstream-checkout", type=Path, help="Original Git checkout, for revision verification and imports")
    parser.add_argument("--dataset-root", type=Path, default=Path("data"))
    parser.add_argument("--region", choices=("belle-ile", "frh04"), default="belle-ile")
    parser.add_argument("--weights", type=Path, default=Path("data/transformer-official.pth"))
    parser.add_argument("--download-official-weights", action="store_true", help="Opt in to downloading published checkpoint")
    parser.add_argument("--offline", action="store_true", help="Refuse any hidden dataset asset downloads")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--max-samples", type=int, default=None, help="Smoke test only; disables historical acceptance threshold")
    parser.add_argument("--output", type=Path, default=Path("results/breizhcrops.json"))
    parser.add_argument("--require-upstream-revision", action="store_true", help="Fail unless checkout matches historical SHA")
    return parser.parse_args()


def run(args):
    if args.batch_size < 1 or any(s < 0 or s > 2**32 - 1 for s in args.seeds):
        raise ValueError("Invalid batch size or seed")
    if args.offline and args.download_official_weights:
        raise ValueError("--offline conflicts with --download-official-weights")
    revision = git_revision(args.upstream_checkout)
    if args.require_upstream_revision and (revision is None or not revision["matches_target"] or revision["dirty"]):
        raise ValueError("Upstream checkout is not the pristine target revision")
    result = {"experiment": "BreizhCrops published Transformer, 2017 L1C", "paper_doi": DOI,
              "upstream_repository": "https://github.com/dl4sits/BreizhCrops", "expected_revision": UPSTREAM_SHA,
              "upstream_checkout": revision, "region": args.region, "level": "L1C", "year": 2017,
              "max_samples": args.max_samples, "seeds": args.seeds, "batch_size": args.batch_size,
              "environment": {"python": platform.python_version(), "torch": torch.__version__, "numpy": np.__version__,
                              "platform": platform.platform(), "cpu_only": True},
              "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "trials": [], "status": "blocked"}
    try:
        if args.upstream_checkout:
            sys.path.insert(0, str(args.upstream_checkout.resolve()))
        if args.download_official_weights:
            download_official_weights(args.weights)
        if not args.weights.is_file():
            raise FileNotFoundError(f"Checkpoint unavailable: {args.weights}; explicitly supply --weights or --download-official-weights")
        result["checkpoint_sha256"] = sha256(args.weights)
        model = load_original_model(args.weights)
        ds_mod = importlib.import_module("breizhcrops.datasets.breizhcrops")
        if args.offline:
            def deny_download(*_a, **_kw):
                raise RuntimeError("Offline mode: required official dataset files are not cached")
            ds_mod.download_file = deny_download
        dataset = ds_mod.BreizhCrops(region=args.region, root=str(args.dataset_root.resolve()),
                                     year=2017, level="L1C", preload_ram=False)
        result["samples_in_dataset"] = len(dataset)
        hashes = {}
        for name in ("h5path", "indexfile", "codesfile", "classmapping"):
            p = Path(getattr(dataset, name))
            if p.is_file():
                hashes[name] = {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}
        result["dataset_files"] = hashes
        for seed in args.seeds:
            result["trials"].append(evaluate(model, dataset, seed=seed,
                                              batch_size=args.batch_size, max_samples=args.max_samples))
        threshold = None if args.max_samples is not None else HISTORICAL_MINIMUM.get((args.region, "L1C"))
        result["historical_minimum"] = threshold
        result["mean_accuracy"] = float(np.mean([r["accuracy"] for r in result["trials"]]))
        result["min_accuracy"] = float(min(r["accuracy"] for r in result["trials"]))
        result["status"] = ("smoke_only" if threshold is None else
                            "pass" if all(r["accuracy"] >= threshold for r in result["trials"]) else "below_historical_threshold")
    except Exception as exc:
        result["status"] = "blocked"
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
    finally:
        atomic_json(args.output, result)
        print(json.dumps({"status": result["status"], "report": str(args.output),
                          "mean_accuracy": result.get("mean_accuracy"), "error": result.get("error")}, indent=2))
    return 0 if result["status"] in ("pass", "smoke_only") else 1


if __name__ == "__main__":
    raise SystemExit(run(build_args()))
