"""Render the 2025 CC-BY manuscript's already-synchronized pressure plots.

The original Figure 5 is used as an EXPERIMENTAL time-registration reference,
never to tune the theoretical capillary pressure model. Source paper is not
committed or republished wholesale. Small manuscript-figure analysis crops
are an ephemeral Actions artifact (7 days).
"""
from pathlib import Path
import urllib.request
import hashlib
import fitz
import json
import io

HERE=Path(__file__).resolve().parent
OUT=HERE/"outputs"
OUT.mkdir(parents=True,exist_ok=True)
SOURCES=[
  "https://pure.tudelft.nl/ws/portalfiles/portal/244708107/1-s2.0-S0021979725012755-main.pdf",
  "https://repository.tudelft.nl/file/File_6f338f49-47aa-4134-b5cf-874dc7654a02",
]
blob=None
errors=[]
for url in SOURCES:
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"AcademicEvidenceAudit/1.0"})
        with urllib.request.urlopen(req,timeout=45) as f:
            raw=f.read(18_000_001)
        if not raw.startswith(b"%PDF") or len(raw)>18_000_000:
            raise ValueError("No valid PDF or over 18MB")
        blob=raw
        source=url
        break
    except Exception as ex:
        errors.append(f"{url}: {type(ex).__name__} {str(ex)[:160]}")
if blob is None:
    (OUT/"paper_fetch_status.json").write_text(json.dumps({"status":"failed","errors":errors},indent=2))
    raise RuntimeError("Unable to acquire CC-BY published manuscript: "+repr(errors))
doc=fitz.open(stream=blob,filetype="pdf")
matches=[]
for j,p in enumerate(doc):
    txt=p.get_text()
    if "Fig. 5." in txt and "pressure" in txt.lower() and "curvature" in txt.lower():
        matches.append(j)
if not matches:
    raise RuntimeError("Published Figure 5 was not found by text search")
# Figure may be cited in neighboring pages: select page with 9-panel figure
index=matches[-1]
pg=doc[index]
pix=pg.get_pixmap(matrix=fitz.Matrix(1.8,1.8),alpha=False)
(png:=OUT/"published_fig5_source_page.png").write_bytes(pix.tobytes("png"))
report={"source":source,"sha256":hashlib.sha256(blob).hexdigest(),
        "pdf_pages":len(doc),"candidate_pages":matches,"rendered_page":index,
        "width_px":pix.width,"height_px":pix.height,
        "status":"source PDF rendered; pressure trace digitization NOT yet attempted"}
(OUT/"publication_figure_provenance.json").write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2),flush=True)
