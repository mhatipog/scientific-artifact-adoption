"""Read selected frames only from a public 329 MB TIFF, via bounded HTTP ranges.

CC BY 4.0: Yan/Müller et al., doi:10.18419/DARUS-4683.
This is exploratory event alignment; camera and pressure time origins differ/are unverified.
"""
import io,json
from pathlib import Path
from urllib.request import Request,urlopen
URL="https://darus.uni-stuttgart.de/api/access/datafile/351438"
SIZE=345104065
OUT=Path(__file__).parent/"outputs"/"empirical"
OUT.mkdir(parents=True,exist_ok=True)
class Remote(io.RawIOBase):
 def __init__(self):
  super().__init__()
  self.pos=0;self.name="120_Q1_water_air_4700-5357_TD.tif"
  self.cache_offset=-1;self.cache=b"";self.calls=0;self.bytes=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,o,whence=0):
  self.pos=int(o if whence==0 else self.pos+o if whence==1 else SIZE+o)
  if not 0<=self.pos<=SIZE:raise ValueError("bad seek")
  return self.pos
 def read(self,n=-1):
  if not 0<=n<=2_000_000:raise RuntimeError("unbounded/oversized read forbidden")
  if self.cache_offset<=self.pos and self.pos+n<=self.cache_offset+len(self.cache):
   q=self.pos-self.cache_offset;self.pos+=n;return self.cache[q:q+n]
  length=min(SIZE-self.pos,max(n,2048) if n<64000 else n)
  if self.calls>=850 or self.bytes+length>24_000_000:
   raise RuntimeError("remote sparse-access budget exceeded")
  if length<=0:return b""
  req=Request(URL,headers={"Range":f"bytes={self.pos}-{self.pos+length-1}",
    "User-Agent":"Open-research-sparse-TIFF/1.0"})
  with urlopen(req,timeout=25) as r:
   if r.status!=206 or not r.headers.get("Content-Range","").startswith(f"bytes {self.pos}-"):
    raise RuntimeError("server does not support sparse range")
   dat=r.read(length+1)
  if len(dat)!=length:raise RuntimeError("wrong byte count")
  self.calls+=1;self.bytes+=len(dat)
  if self.calls%100==0:print("REMOTE RANGE PROGRESS",self.calls,self.bytes,flush=True)
  self.cache_offset=self.pos;self.cache=dat;self.pos+=n;return dat[:n]
 def readinto(self,b):
  view=memoryview(b).cast('B');v=self.read(view.nbytes);view[:len(v)]=v;return len(v)
def main():
 import numpy as np,tifffile
 from PIL import Image
 r=Remote()
 report={"status":"started","file_id":351438,"source_bytes":SIZE,
  "requested_indices":[0,1,5,10,20,30,50,75,100,150,200,250,320,340,360,657],
  "warning":"video and pressure logger time origins UNVERIFIED"}
 try:
  with tifffile.TiffFile(r) as tf:
   report["pages"]=len(tf.pages)
   report["first_page_shape"]=list(tf.pages[0].shape)
   for index in report["requested_indices"]:
    if index>=len(tf.pages):continue
    f=tf.pages[index].asarray()
    if f.ndim>2:f=f[...,0]
    p1,p99=np.percentile(f,[1,99]);p99=max(p99,p1+1)
    g=np.clip((f.astype(float)-p1)/(p99-p1)*255,0,255).astype("uint8")
    p=OUT/"yan_frames";p.mkdir(exist_ok=True)
    name=f"120_Q1_{index:03d}_frame_{4700+index}.png"
    Image.fromarray(g).resize((512,512)).save(p/name,optimize=True)
    inside=f[115:290,120:280].astype(float)
    background=f[365:435,200:330].astype(float)
    report.setdefault("frames",[]).append({"index":index,"original_frame":4700+index,
      "hypothetical_camera_elapsed_s":(4700+index)/30,
      "shape":list(f.shape),"dtype":str(f.dtype),"png":name,
      "pore_median_raw":float(np.median(inside)),
      "background_median_raw":float(np.median(background)),
      "pore_to_background_ratio":float(np.median(inside)/(np.median(background)+1e-12))})
  report["status"]="success_sparse_frames"
 except Exception as e:
  report["status"]="blocked_sparse_access";report["error"]=repr(e)[:600]
 report["requests"]=r.calls;report["bytes_received"]=r.bytes
 (OUT/"yan_sparse_frame_report.json").write_text(json.dumps(report,indent=2))
 print("SPARSE IMAGE REPORT",json.dumps(report),flush=True)
if __name__=="__main__":main()
