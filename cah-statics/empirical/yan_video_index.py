"""Video archive metadata; no TIFF download or invented camera synchronization."""
import re,json
from pathlib import Path
from urllib.request import urlopen,Request
from urllib.parse import quote
url="https://darus.uni-stuttgart.de/api/datasets/:persistentId/?persistentId="+quote("doi:10.18419/DARUS-4683",safe="")
with urlopen(Request(url,headers={"User-Agent":"Independent-CAH-Reproducibility/1.0"}),timeout=35) as r:
    metadata=json.load(r)
files=metadata["data"]["latestVersion"]["files"]
out=[]
for row in files:
    o=row["dataFile"]
    name=o["filename"]
    if not name.endswith(".tif") or "_Q1_water_air_" not in name:continue
    m=re.search(r'^(\d+)_Q1_water_air_(\d+)-(\d+).*\.tif$',name,re.I)
    if not m:continue
    angle,lo,hi=map(int,m.groups())
    out.append(dict(angle_deg=angle,filename=name,file_id=o.get("id"),
                    frames=hi-lo+1,relative_frame_window_s=[round(lo/30,3),round(hi/30,3)],
                    video_reference="frame numbering assumes 30fps and shared camera origin; pressure logger origin UNCONFIRMED"))
out.sort(key=lambda x:x["angle_deg"])
path=Path(__file__).parent/"outputs"/"empirical"
path.mkdir(parents=True,exist_ok=True)
(path/"yan_video_metadata_no_download.json").write_text(json.dumps(out,indent=2))
print("Archive low-flow water-air video index: ",json.dumps(out))
print("CAUTION: video frame numbers are not proven synchronized to pressure's t=0.")
