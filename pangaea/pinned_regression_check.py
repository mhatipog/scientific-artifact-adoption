"""Check PANGAEA's ORIGINAL regression evaluator under alternative execution.

This script fetches NOTHING: point --repo at an independently checked-out
PANGAEA commit. It extracts the original RegEvaluator.evaluate AST unchanged,
and extracts a minimal correction made only to its three MSE accumulation
operations. In both cases the actual model/evaluator code executes with real
PyTorch distributed collectives. No BioMassters checkpoint is used.
"""
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

import torch
import torch.nn.functional as F

PINNED_COMMIT = "ca19b837a16797cd809225225b671473c7c93768"
PINNED_BLOB = "345e7b33a471134e2d1a5bd3f1d28c9d5776b73c"
SOURCE = Path("pangaea/engine/evaluator.py")
CHANGES = (
    (
        "        mse = torch.zeros(1, device=self.device)\n",
        "        error_sum_and_count = torch.zeros(2, device=self.device, dtype=torch.float64)\n"
    ),
    (
        "            mse += F.mse_loss(logits, target)\n",
        "            error_sum_and_count[0] += F.mse_loss(\n"
        "                logits.double(), target.double(), reduction='sum'\n"
        "            )\n"
        "            error_sum_and_count[1] += target.numel()\n"
    ),
    (
        "        torch.distributed.all_reduce(mse, op=torch.distributed.ReduceOp.SUM)\n"
        "        mse = mse / len(self.val_loader)\n",
        "        torch.distributed.all_reduce(error_sum_and_count, op=torch.distributed.ReduceOp.SUM)\n"
        "        if error_sum_and_count[1].item() == 0:\n"
        "            raise ValueError('No evaluated pixels')\n"
        "        mse = error_sum_and_count[0] / error_sum_and_count[1]\n"
    ),
)

def correct(source):
    for old, new in CHANGES:
        if source.count(old) != 1:
            raise AssertionError("Pinned method changed; inspect source before patch: "+repr(old))
        source = source.replace(old, new, 1)
    return source

def extract(source, filename):
    tree = ast.parse(source, filename=filename)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef)
               and n.name == "RegEvaluator")
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                  and n.name == "evaluate")
    cloned = copy.deepcopy(method)
    audit_cls = ast.ClassDef(name="OriginalMethod", bases=[], keywords=[],
                             body=[cloned], decorator_list=[])
    module = ast.fix_missing_locations(ast.Module(body=[audit_cls], type_ignores=[]))
    # No mocked arithmetic: only the unused PANGAEA dependency graph is
    # excluded from this execution of the *unmodified evaluator method*.
    namespace = dict(torch=torch, F=F, time=time, os=os,
                     tqdm=lambda sequence, **kw: sequence)
    exec(compile(module, filename, "exec"), namespace)
    return namespace["OriginalMethod"]

class FixedPredictionModel:
    def eval(self):
        return self
    def __call__(self, image, output_shape=None):
        return image["x"]

def score(evaluator_type, values, batch_size):
    batches = []
    for i in range(0,len(values),batch_size):
        batch = torch.tensor(values[i:i+batch_size], dtype=torch.float32)
        batches.append(dict(image=dict(x=batch[:,None,None,None]),
                            target=torch.zeros(batch.shape[0],1,1)))
    obj = evaluator_type()
    obj.device = torch.device("cpu")
    obj.split = "test"
    obj.inference_mode = "whole"
    obj.val_loader = batches
    obj.log_metrics = lambda metrics: None
    metrics, _ = obj.evaluate(FixedPredictionModel())
    return metrics["MSE"]

def provenance(repo):
    repo = repo.resolve()
    head = subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],
                                   text=True).strip()
    sha = subprocess.check_output(["git","-C",str(repo),"hash-object",str(repo/SOURCE)],
                                  text=True).strip()
    if head != PINNED_COMMIT or sha != PINNED_BLOB:
        raise AssertionError(f"Unexpected PANGAEA revision: {head} {sha}")
    source = (repo/SOURCE).read_text()
    ast.parse(source)
    print("PROVENANCE",json.dumps(dict(commit=head,blob=sha,
               torch=torch.__version__,file=str(repo/SOURCE))),flush=True)
    return source

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path,required=True)
    ap.add_argument("--mode",choices=("single","distributed"),required=True)
    ap.add_argument("--results",type=Path,default=Path("pangaea-regression-results.json"))
    args=ap.parse_args()
    source=provenance(args.repo)
    original=extract(source,str(args.repo/SOURCE))
    patched=extract(correct(source),str(args.repo/SOURCE))
    if args.mode=="single":
        with tempfile.TemporaryDirectory() as d:
            torch.distributed.init_process_group("gloo",
                        init_method="file://"+str(Path(d)/"init"),rank=0,world_size=1)
            try:
                result=[]
                for batch in (1,3,8,9):
                    a=[0.0]*8+[3.0]
                    b=[2.0**0.5]*9
                    item=dict(batch_size=batch,
                              original_A=score(original,a,batch),
                              original_B=score(original,b,batch),
                              corrected_A=score(patched,a,batch),
                              corrected_B=score(patched,b,batch))
                    result.append(item)
                assert all(abs(x["corrected_A"]-1)<1e-5 and
                           abs(x["corrected_B"]-2)<1e-5 for x in result)
                assert result[0]["original_A"]<result[0]["original_B"]
                assert result[2]["original_A"]>result[2]["original_B"]
                print("SINGLE_PROCESS_INVARIANCE",json.dumps(result),flush=True)
                args.results.write_text(json.dumps(result,indent=2)+"\n")
            finally:
                torch.distributed.destroy_process_group()
    else:
        torch.distributed.init_process_group("gloo",init_method="env://")
        try:
            rank=torch.distributed.get_rank()
            if torch.distributed.get_world_size()!=2: raise AssertionError("Need two ranks")
            preds=[float(n)**0.5 for n in range(1+4*rank,5+4*rank)]
            before=score(original,preds,4)
            after=score(patched,preds,4)
            if rank==0:
                output=dict(world_size=2,old_mse=before,corrected_mse=after,
                            mathematically_correct_mse=4.5)
                assert abs(before-9)<1e-5,output
                assert abs(after-4.5)<1e-5,output
                print("TWO_RANK_INVARIANCE",json.dumps(output),flush=True)
                args.results.write_text(json.dumps(output,indent=2)+"\n")
        finally:
            torch.distributed.destroy_process_group()

if __name__=="__main__":
    main()
