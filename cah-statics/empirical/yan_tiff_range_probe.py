"""Safe HTTP byte-range probe for the 120° Q1 329MB TIFF.

Never fetches the full file. Stops after 16 KiB regardless of response.
Does not interpret a partial TIFF as an experimental frame.
"""
import json
from pathlib import Path
from urllib.request import Request,urlopen
url="https://darus.uni-stuttgart.de/api/access/datafile/351438"
out=Path(__file__).parent/"outputs"/"empirical"
out.mkdir(parents=True,exist_ok=True)
results={"source_url":url,"purpose":"can remote TIFF be randomly accessed without huge data transfer?",
         "full_video_downloaded":False}
try:
    req=Request(url,headers={"Range":"bytes=0-16383",
        "User-Agent":"Research-Sparse-Frame-Access/1.0"})
    with urlopen(req,timeout=35) as f:
        status=f.status
        content_range=f.headers.get("Content-Range")
        accept_ranges=f.headers.get("Accept-Ranges")
        length=f.headers.get("Content-Length")
        sample=f.read(16384)
    results.update(http_status=status,content_range=content_range,
                   accept_ranges=accept_ranges,reported_content_length=length,
                   sample_bytes_read=len(sample),header_signature=sample[:4].hex(),
                   safe_random_access=(status==206 and content_range is not None))
except Exception as e:results["error"]=str(e)[:250]
(out/"yan_120_tiff_range_probe.json").write_text(json.dumps(results,indent=2))
print("TIFF RANGE PROBE",json.dumps(results))
