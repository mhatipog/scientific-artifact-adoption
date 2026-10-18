"""Retrieve ONE public microscope TIFF and produce a *derived* contact sheet.

Data: Yan/Müller et al. 2025, DaRUS DOI 10.18419/DARUS-4683 (CC BY 4.0).
This does not label pore-filling times or claim pressure/model agreement.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen
from urllib.parse import quote
import json,re,hashlib,sys
import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE="https://darus.uni-stuttgart.de"
PID="doi:10.18419/DARUS-4683"
HERE=Path(__file__).resolve().parent
OUT=HERE/"outputs"
OUT.mkdir(exist_ok=True)
LIMIT=430_000_000
HEADERS={"User-Agent":"IndependentCAHScientificAudit/1.0","Accept":"*/*"}

def open_url(url,timeout=90):
    return urlopen(Request(url,headers=HEADERS),timeout=timeout)

def scan_tiff(path,chosen):
    with tifffile.TiffFile(str(path)) as t:
        ser=t.series[0];pages=t.pages
        info={"series_shape":list(ser.shape),"axes":ser.axes,"page_count":len(pages),
              "first_page_shape":list(pages[0].shape),"imagej":bool(t.is_imagej),
              "first_page_tags":{key:str(pages[0].tags[key].value)[:300]
                 for key in ("DateTime","ImageDescription","XResolution","YResolution","ResolutionUnit")
                 if key in pages[0].tags}}
        if len(pages)<2:
            raise RuntimeError("TIFF did not expose multiple independently readable pages")
        indices=np.unique(np.rint(np.linspace(0,len(pages)-1,12)).astype(int))
        frames=[]
        for ix in indices:
            a=pages[int(ix)].asarray()
            if a.ndim==3 and a.shape[-1] in (3,4):a=a[...,:3].mean(axis=-1)
            if a.ndim!=2 or max(a.shape)>4096:raise RuntimeError(f"unsupported frame dimensions {a.shape}")
            frames.append(np.asarray(a,dtype=np.float32))
    lo,hi=np.percentile(np.stack(frames),[1,99.5])
    fig,axes=plt.subplots(3,4,figsize=(13,9),layout="constrained")
    for ax,ix,f in zip(axes.flat,indices,frames):
        ax.imshow(f,cmap="gray",vmin=lo,vmax=hi,interpolation="nearest")
        ax.set(title=f"Frame {ix} / {info['page_count']}")
        ax.axis("off")
    fig.suptitle("Original microscopy: "+chosen["filename"]+"\n"
                 "Yan/Müller et al. (2025), DOI 10.18419/DARUS-4683, CC BY 4.0\n"
                 "Original frames, no pressure-event alignment or model validation",fontsize=13)
    fig.savefig(OUT/"microscopy_contact_sheet.png",dpi=120)
    plt.close(fig)
    info["sample_indices"]=indices.tolist()
    info["sample_median_change_from_first"]=[round(float(np.median(abs(f-frames[0]))),4) for f in frames]
    info["sample_mean_intensity"]=[round(float(np.mean(f)),4) for f in frames]
    return info

def main():
    report={"dataset":PID,"status":"not started"}
    try:
        url=BASE+"/api/datasets/:persistentId/?persistentId="+quote(PID,safe="")
        with open_url(url,45) as r: meta=json.load(r)
        files=[a["dataFile"] for a in meta["data"]["latestVersion"]["files"]]
        options=[{"filename":f["filename"],"id":f["id"],"bytes":f["filesize"],
                  "expected_md5":f.get("md5",f.get("checksum",{}).get("value"))}
                 for f in files if re.fullmatch(r"(60|95|120)_Q1_water_air_[0-9-]+(?:_TD)?\.tif",f["filename"],re.I)]
        options.sort(key=lambda a:a["bytes"])
        report["available_q1_water_air_tiffs"]=options
        print("TIFF_OPTIONS",json.dumps(options),flush=True)
        feasible=[q for q in options if q["bytes"]<LIMIT]
        if not feasible:raise RuntimeError("No Q1 water-air TIFF under safe 430 MB download cap")
        chosen=next((q for q in feasible if q["filename"].startswith("95_Q1_")),feasible[0])
        report["chosen"]=chosen
        with TemporaryDirectory(prefix="cah-yan-video-") as td:
            path=Path(td)/chosen["filename"]
            md5=hashlib.md5();total=0
            with open_url(BASE+"/api/access/datafile/"+str(chosen["id"]),90) as response,open(path,"wb") as target:
                if "html" in response.headers.get("Content-Type","").lower():
                    raise RuntimeError("HTML returned, not a TIFF")
                while True:
                    data=response.read(2_000_000)
                    if not data:break
                    total+=len(data)
                    if total>LIMIT:raise RuntimeError("Exceeded size cap")
                    target.write(data);md5.update(data)
            if chosen["expected_md5"] and md5.hexdigest()!=chosen["expected_md5"]:
                raise RuntimeError("Dataset MD5 mismatch")
            report["download"]={"bytes":total,"verified_md5":md5.hexdigest()}
            report["microscopy"]=scan_tiff(path,chosen)
            from track_video import measure
            report["optical_filling"]=measure(str(path),chosen)
        report["status"]="original TIFF inspected; derived contact sheet made; synchronization outstanding"
    except Exception as error:
        report["status"]="blocked";report["error"]=repr(error)
        print("BLOCKED",repr(error),flush=True)
    (OUT/"video_metadata.json").write_text(json.dumps(report,indent=2))
    print("RESULT",json.dumps(report,indent=2)[:7000],flush=True)
    if report["status"]=="blocked":sys.exit(1)

if __name__=="__main__":main()
