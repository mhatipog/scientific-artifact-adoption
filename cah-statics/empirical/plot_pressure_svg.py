"""Self-contained, dependency-free SVG of DaRUS measured pressure traces.

Runs only *after* audit_darus.py. Raw measurements do not enter git.
SVG is a graphical reading aid, not an event-labelled model validation.
"""
from pathlib import Path
import json, html

OUT=Path(__file__).parent/"outputs"/"empirical"
DATA=json.loads((OUT/"darus_lowflow_pressure_audit.json").read_text())
C=["#027d98","#d2673d","#7651a3"]
w,h=1240,920
panels=[]
S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
   '<rect width="1240" height="920" fill="#fbfdff"/>',
   '<text x="60" y="55" font-family="Arial" font-weight="700" font-size="27" fill="#183047">Independent measured pressure recordings | Yan–Müller 2025</text>',
   '<text x="60" y="82" font-family="Arial" font-size="15" fill="#55667b">Original sensor files, kPa; pressure-event windows are not aligned to video or predicted by the model.</text>']
for i,(f,col) in enumerate(zip(DATA["files"],C)):
    vals=f.get("decimated_trace",[])
    if len(vals)<3:continue
    tmax=max(a[0] for a in vals)
    ymin=min(a[1] for a in vals);ymax=max(a[1] for a in vals)
    pad=max(.15,(ymax-ymin)*.08);ymin-=pad;ymax+=pad
    x0,y0=95,150+235*i
    tw,th=1070,152
    proj=lambda t,p:(x0+tw*t/tmax,y0+th*(ymax-p)/(ymax-ymin))
    d=" ".join(("M" if j==0 else "L")+f"{proj(t,p)[0]:.2f},{proj(t,p)[1]:.2f}" for j,(t,p) in enumerate(vals))
    for j in range(5):
        Y=y0+th*j/4
        yy=ymax-(ymax-ymin)*j/4
        S.append(f'<line x1="{x0}" x2="{x0+tw}" y1="{Y:.2f}" y2="{Y:.2f}" stroke="#dce3ed" stroke-width="1"/>')
        S.append(f'<text x="86" y="{Y+4:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#63778b">{yy:.2f}</text>')
    S.append(f'<path d="{d}" stroke="{col}" fill="none" stroke-width="2.15"/>')
    ev=f.get("largest_five_second_decline_native",{})
    if ev:
        xstart=proj(ev["t_start_s"],0)[0];xend=proj(ev["t_end_s"],0)[0]
        S.append(f'<rect x="{xstart:.1f}" y="{y0}" width="{max(2,xend-xstart):.1f}" height="{th}" fill="#eec562" opacity=".37"/>')
    theta=f["name"].split("_")[0]
    S.append(f'<text x="{x0}" y="{y0-27}" font-family="Arial" font-weight="700" font-size="17" fill="#183047">{theta}° · {f.get("valid_samples",0):,} samples · {f.get("elapsed_seconds",0)/60:.1f} min</text>')
    S.append(f'<text x="{x0+tw}" text-anchor="end" y="{y0-27}" font-family="Arial" font-size="13" fill="#64748b">Largest 5-second decline {ev.get("delta",float("nan")):.3f} kPa (unverified event)</text>')
    for j in range(5):
        xx=x0+tw*j/4;tt=tmax*j/4
        S.append(f'<text x="{xx:.1f}" y="{y0+th+22}" font-family="Arial" font-size="12" text-anchor="middle" fill="#64748b">{tt:.0f}s</text>')
    S.append(f'<text x="{x0+tw/2}" y="{y0+th+42}" font-family="Arial" font-size="12" text-anchor="middle" fill="#64748b">Time since start of this pressure recording</text>')
S.append('<text x="60" y="884" font-family="Arial" font-size="13" fill="#55667b">Source: Müller, Yan, Raoof &amp; Weigand (2025), DaRUS doi:10.18419/DARUS-4683 (CC BY 4.0). </text>')
S.append('<text x="60" y="905" font-family="Arial" font-size="13" fill="#9c4c41">CAUTION: displayed declines are signal diagnostics, NOT individually identified capillary jumps or quantitative model agreement.</text>')
S.append("</svg>")
(OUT/"measured_pressure_traces.svg").write_text("\n".join(S))
print("PASS SVG from authentic 60/95/120-degree decimated timestamped pressure data; no jump fitting")
