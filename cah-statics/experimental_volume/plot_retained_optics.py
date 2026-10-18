"""Deterministic, offline SVG from frozen 95-degree seven-pixel optical data.

This plots dimensionless optical change, NOT liquid thickness/volume.
Original data: Yan et al. (2025), doi:10.18419/DARUS-4683, CC BY 4.0.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from xml.sax.saxutils import escape

DATA=Path(__file__).resolve().parent/"derived"/"retained_optical_95_4822_4824_7px.csv"
OUT=Path(__file__).resolve().parent/"figures"/"retained_optical_progress_7px.svg"


def svg_string(rows):
    if len(rows)!=5 or any(int(r["retained_boundary_exclusion_px"])!=7 for r in rows):
        raise ValueError("Expected exact five seven-pixel boundary-exclusion rows")
    width,height=890,490
    left,right,top,bottom=88,850,95,384
    ymax=.22
    def Y(v):
        return bottom-(float(v)/ymax)*(bottom-top)
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Normalized optical change among retained-wet and stable-off pixels">',
       '<rect width="890" height="490" fill="white"/>',
       '<text x="88" y="37" font-family="Arial,sans-serif" font-size="21" font-weight="700" fill="#162534">Original 95° frames 4822–4824: pixels retained wet still change optically</text>',
       '<text x="88" y="61" font-family="Arial,sans-serif" font-size="14" fill="#465867">Seven-pixel boundary exclusion (~17.4 μm); five frozen classifiers; NOT liquid depth</text>']
    for v in (0,.05,.10,.15,.20):
        y=Y(v)
        s.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#e2e9ed"/>')
        s.append(f'<text x="{left-10}" y="{y+5:.2f}" font-family="Arial,sans-serif" font-size="13" text-anchor="end" fill="#465867">{v:.2f}</text>')
    for i,r in enumerate(rows):
        x=152+i*149
        retain=float(r["retained_delta_median"])
        off=float(r["stable_off_delta_median"])
        for dx,v,color in ((0,retain,"#2767a3"),(33,off,"#b87438")):
            y=Y(v)
            s.append(f'<rect x="{x+dx}" y="{y:.2f}" width="25" height="{bottom-y:.2f}" fill="{color}"/>')
        label=escape(str(r["count_threshold"])+"/"+str(r["fraction_threshold"]))
        s.append(f'<text x="{x+29}" y="408" font-family="Arial,sans-serif" font-size="14" text-anchor="middle" fill="#273844">{label}</text>')
        s.append(f'<text x="{x+29}" y="427" font-family="Arial,sans-serif" font-size="11" text-anchor="middle" fill="#61707e">n={int(r["retained_wet_count"]):,}</text>')
    s.extend([
        '<text x="88" y="458" font-family="Arial,sans-serif" font-size="12" fill="#61707e">Classifier: reference-to-final contrast counts / fractional wet threshold; retained counts n are spatially correlated</text>',
        '<rect x="560" y="68" width="14" height="14" fill="#2767a3"/>',
        '<text x="581" y="80" font-family="Arial,sans-serif" font-size="13" fill="#273844">Retained-wet median Δr</text>',
        '<rect x="560" y="86" width="14" height="14" fill="#b87438"/>',
        '<text x="581" y="98" font-family="Arial,sans-serif" font-size="13" fill="#273844">Stable-off median Δr</text>',
        '</svg>'])
    return "\n".join(s)+"\n"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stdout",action="store_true")
    args=parser.parse_args()
    with DATA.open(newline="") as f:
        rows=list(csv.DictReader(f))
    figure=svg_string(rows)
    if args.stdout:
        print(figure,end="")
    else:
        OUT.parent.mkdir(exist_ok=True,parents=True)
        OUT.write_text(figure,encoding="utf-8")
        print(OUT)


if __name__=="__main__":
    main()
