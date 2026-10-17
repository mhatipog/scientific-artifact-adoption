"""Cross-check actual PANGAEA segmentation and multilabel evaluator metrics.

The evaluator methods are AST-extracted unchanged from pinned original code.
Only unrelated imports and dataset/model infrastructure are stubbed.
"""
import argparse
import ast
import copy
import json
import logging
import os
from pathlib import Path
import time
import torch
import sklearn
import sklearn.metrics
from pinned_regression_check import provenance, SOURCE

def extract(source, clsname, method_name):
    original = ast.parse(source)
    cls = next(x for x in original.body if isinstance(x,ast.ClassDef)
               and x.name==clsname)
    method = copy.deepcopy(next(x for x in cls.body if isinstance(x,ast.FunctionDef)
                                and x.name==method_name))
    standin=ast.ClassDef(name="ActualMethod",bases=[],keywords=[],
                        body=[method],decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[standin],type_ignores=[]))
    namespace=dict(torch=torch,sklearn=sklearn,time=time,os=os,
                   Path=Path,KNNClassifier=object,DataLoader=object,
                   tqdm=lambda seq,**kwargs: seq)
    exec(compile(module,str(SOURCE),"exec"),namespace)
    return namespace["ActualMethod"]

def segmentation(source):
    cls=extract(source,"SegEvaluator","compute_metrics")
    cm=torch.tensor([[40,10],[20,30]],dtype=torch.float32)
    obj=cls()
    obj.valid_class_indices=[0,1]
    res=obj.compute_metrics(cm)
    true_labels=[0]*60+[1]*40
    predicted_labels=[0]*40+[1]*20+[0]*10+[1]*30
    check=torch.bincount(torch.tensor(predicted_labels)*2+
                         torch.tensor(true_labels),minlength=4).view(2,2)
    assert torch.all(check==cm.long())
    expected=sklearn.metrics.precision_recall_fscore_support(
        true_labels,predicted_labels,labels=[0,1],zero_division=0)
    assert all(abs(res["Precision"][i]/100-expected[1][i])<1e-6 for i in range(2))
    assert all(abs(res["Recall"][i]/100-expected[0][i])<1e-6 for i in range(2))
    assert any(abs(res["Precision"][i]/100-expected[0][i])>.05 for i in range(2))
    output=dict(confusion_matrix_rows_pred_cols_truth=cm.int().tolist(),
                upstream_precision=res["Precision"],upstream_recall=res["Recall"],
                sklearn_precision=expected[0].tolist(),
                sklearn_recall=expected[1].tolist())
    print("ACTUAL_SEGMENTATION_METRICS",json.dumps(output,sort_keys=True),flush=True)
    return output

class FrozenKnnOutput:
    def __init__(self):
        self.topk=[1]
        self._bank=torch.tensor([1])
        self._bank_labels=torch.tensor([1])
        self.num_classes=3
    def eval(self):
        return self
    def classify(self,image):
        return image["pred"]

def knn_score(cls,targets,preds,bs):
    instance=cls()
    instance.val_loader=[
       dict(image=dict(pred=torch.tensor(preds[i:i+bs],dtype=torch.int64)),
            target=torch.tensor(targets[i:i+bs],dtype=torch.int64))
       for i in range(0,len(targets),bs)]
    instance.device=torch.device("cpu")
    instance.logger=logging.getLogger("audit")
    instance.split="test"
    instance.multi_label=True
    instance.use_wandb=False
    result,_=instance.evaluate(FrozenKnnOutput(),train_loader=None,model_name="knn")
    return float(result["f1"])

def multilabel(source):
    cls=extract(source,"KNNClassificationEvaluator","evaluate")
    actual=[[1,1,1],[0,0,1],[1,0,0],[0,0,0]]
    pred=[[1,1,1],[0,0,1],[0,0,0],[1,1,1]]
    f1_1=knn_score(cls,actual,pred,1)
    f1_4=knn_score(cls,actual,pred,4)
    reference=float(sklearn.metrics.f1_score(actual,pred,average="micro",zero_division=0))
    assert abs(reference-2/3)<1e-9
    assert abs(f1_1-0.5)<1e-9
    assert abs(f1_4-reference)<1e-9
    output=dict(batch_1=f1_1,batch_4=f1_4,global_micro_f1=reference)
    print("ACTUAL_MULTILABEL_KNN",json.dumps(output,sort_keys=True),flush=True)
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--repo",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    source=provenance(args.repo)
    obj=dict(segmentation=segmentation(source),multi_label_knn=multilabel(source))
    args.out.write_text(json.dumps(obj,indent=2)+"\n")
    print("SOURCE_LEVEL_OTHER_METRICS_PASSED",flush=True)

if __name__=="__main__":
    main()
