"""Fetch ONLY three CC-BY pressure traces from DaRUS and audit the raw structure.

Never silently treat pressure data as a prediction of a two-meniscus model.
Network access is only during manually triggered / path-scoped PR Action.
No raw pressure files committed to git or uploaded as artifacts.
"""
from __future__ import annotations
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import quote
from collections import Counter
import csv, json, hashlib, re, statistics, traceback

BASE="https://darus.uni-stuttgart.de"
DATASET="doi:10.18419/DARUS-4683"
OUT=Path(__file__).parent/"outputs"/"empirical"
OUT.mkdir(parents=True,exist_ok=True)
NAMES=[f"{angle}_Q1_pressure_water_air.txt" for angle in (60,95,120)]

def get(url,limit=16_000_000,timeout=35):
    req=Request(url,headers={"User-Agent":"AcademicReproducibility/1.0 (+https://github.com/mhatipog/scientific-artifact-adoption)","Accept":"application/json,text/plain,*/*"})
    with urlopen(req,timeout=timeout) as r:
        status=r.status
        content=r.read(limit+1)
        typ=r.headers.get("Content-Type","")
    if len(content)>limit:
        raise RuntimeError("Unexpectedly large dataset; stop.")
    if status!=200:raise RuntimeError(f"HTTP {status}")
    if content[:100].lower().find(b"<html")>=0:raise RuntimeError("HTML returned rather than dataset")
    return content,typ

def main():
    report={"dataset_doi":DATASET,"source_url":BASE+"/dataset.xhtml?persistentId="+quote(DATASET),
            "license":"CC BY 4.0 as stated by original deposit","requested":NAMES,
            "status":"not analyzed yet","files":[]}
    try:
        url=BASE+"/api/datasets/:persistentId/?persistentId="+quote(DATASET,safe="")
        data,typ=get(url,limit=3_000_000)
        metadata=json.loads(data)
        files=metadata["data"]["latestVersion"]["files"]
        byname={v["dataFile"]["filename"]:v["dataFile"] for v in files}
        report["dataset_file_count"]=len(byname)
        for name in NAMES:
            v=byname.get(name)
            if v is None:
                report["files"].append({"name":name,"status":"not found in dataset metadata"})
                continue
            record={"name":name,"file_id":v["id"],"metadata_md5":v.get("md5",v.get("checksum",{}).get("value")),
                    "source_download":BASE+"/api/access/datafile/"+str(v["id"])}
            try:
                raw,mime=get(record["source_download"],limit=8_000_000)
                record["bytes"]=len(raw);record["md5"]=hashlib.md5(raw).hexdigest()
                record["mime"]=mime
                record["sample_first_8_lines"]=raw.decode("utf-8-sig",errors="replace").splitlines()[:8]
                lines=raw.decode("utf-8-sig",errors="replace").splitlines()
                record["total_lines"]=len(lines)
                numeric=[]
                counts=Counter()
                # Keep only a small bounded parsed preview, not the original dataset.
                for line in lines:
                    l=line.strip()
                    if not l or l.startswith(("#","%","//")):continue
                    vals=l.replace(";", " ").replace("\t"," ").split()
                    try:
                        floats=[float(v.replace(",",".")) for v in vals]
                    except ValueError:continue
                    counts[len(floats)]+=1
                    if len(floats)>=2:numeric.append(floats[:4])
                record["numeric_column_widths"]=dict(counts)
                if numeric:
                    columns=list(zip(*[v[:min(map(len,numeric))] for v in numeric]))
                    record["numeric_rows"]=len(numeric)
                    record["first_numeric_rows"]=numeric[:4]
                    record["column_minmax"]=[{"min":min(c),"max":max(c)} for c in columns]
                record["status"]="downloaded and parsed format; NOT model-validated"
            except Exception as ex:
                record["status"]="download failed"
                record["error"]=str(ex)[:400]
            report["files"].append(record)
        report["status"]="metadata obtained; see individual file statuses"
    except Exception as e:
        report["status"]="blocked: dataset metadata unavailable"
        report["error"]=str(e)[:500]
    (OUT/"darus_lowflow_pressure_audit.json").write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps({"status":report["status"],"files":[
       {k:x.get(k) for k in ("name","status","numeric_rows","error","file_id")}
       for x in report["files"]]},indent=2))
    print("No experimental quantitative fitting or agreement asserted.")

if __name__=="__main__":
    main()
