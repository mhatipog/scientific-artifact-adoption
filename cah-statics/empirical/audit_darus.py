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
                print("RAW_HEADER",name,[repr(v[:200]) for v in lines[:6]])
                print("RAW_DATA_LINE",name,[repr(v[:200]) for v in lines[24:75]])
                print("RAW_TRAILER",name,[repr(v[:120]) for v in lines[-5:]])
                from datetime import datetime
                from math import isfinite
                marker=next((i for i,line in enumerate(lines) if line.strip()=="[Data]"),-1)
                if marker<0:raise ValueError("Missing pressure logger [Data] block")
                metadata={}
                for line in lines[:marker]:
                    pair=line.split("\t",1)
                    if len(pair)==2: metadata[pair[0].strip()]=pair[1].strip()
                record["original_instrument_metadata"]={
                    key:value for key,value in metadata.items()
                    if any(t in key.lower() for t in ("unit","phenomenon","model","calibration","sample"))}
                t0=None
                pressures=[]
                times=[]
                invalid=0
                for line in lines[marker+1:]:
                    row=line.split("\t")
                    if len(row)<4:continue
                    if row[3].strip()!="Valid":
                        invalid+=1
                        continue
                    try:
                        time=datetime.fromisoformat(row[0]+"T"+row[1])
                        value=float(row[2])
                        if not isfinite(value):continue
                    except ValueError:continue
                    if t0 is None:t0=time
                    times.append((time-t0).total_seconds())
                    pressures.append(value)
                if len(pressures)<100:raise ValueError("Fewer than 100 timestamped pressure samples")
                from statistics import median
                sampling=median([times[i+1]-times[i] for i in range(min(3000,len(times)-1))])
                record["native_unit_warning"]="Read AlgorithmUnit in metadata, independently verify exported column units before Pa conversion"
                record["valid_samples"]=len(pressures)
                record["invalid_samples"]=invalid
                record["elapsed_seconds"]=round(times[-1],6)
                record["median_sampling_interval_seconds"]=round(sampling,8)
                sorted_p=sorted(pressures)
                record["pressure_native_min"]=min(pressures)
                record["pressure_native_max"]=max(pressures)
                record["pressure_native_q05"]=sorted_p[int(.05*(len(sorted_p)-1))]
                record["pressure_native_q50"]=sorted_p[len(sorted_p)//2]
                record["pressure_native_q95"]=sorted_p[int(.95*(len(sorted_p)-1))]
                record["first_valid_samples"]=list(zip(times[:4],pressures[:4]))
                # One-second blocks suppress high-frequency sensor jitter;
                # the strongest 5s pressure decrease is *not* by itself a pore event.
                from statistics import mean
                block=max(1,round(1/sampling))
                bt=[]; bp=[]
                for j in range(0,len(pressures)-block+1,block):
                    bt.append(mean(times[j:j+block]))
                    bp.append(mean(pressures[j:j+block]))
                lag=5
                if len(bp)>lag:
                    event=min(range(lag,len(bp)),key=lambda j: bp[j]-bp[j-lag])
                    record["largest_five_second_decline_native"]={
                        "delta":round(bp[event]-bp[event-lag],6),
                        "t_start_s":round(bt[event-lag],3),
                        "t_end_s":round(bt[event],3),
                        "classification":"UNALIGNED diagnostic only; not independently a pore-filling jump"}
                record["native_trace_one_second_blocks"]=len(bp)
                # Store at most 900 samples per trace, not original raw observations.
                stride=max(1,len(bp)//900)
                record["decimated_trace"]=[[round(bt[j],2),round(bp[j],5)]
                                            for j in range(0,len(bp),stride)]
                record["status"]="raw pressure series parsed; event alignment and theoretical pressure predictions outstanding"

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
       {k:x.get(k) for k in ("name","status","valid_samples","elapsed_seconds","pressure_native_min","pressure_native_max","largest_five_second_decline_native","error","file_id")}
       for x in report["files"]]},indent=2))
    print("No experimental quantitative fitting or agreement asserted.")

if __name__=="__main__":
    main()
