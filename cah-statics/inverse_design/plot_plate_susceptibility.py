"""Deterministic dependency-light SVG figure for Round 17 FEM convergence."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from plate_susceptibility import report


def make_svg() -> str:
    data=report()
    width,height=760,465
    left,top,right,bottom=98,86,680,362
    panels=[]
    panels.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    panels.append('<rect width="760" height="465" fill="#ffffff"/>')
    panels.append('<text x="66" y="34" font-size="20" font-family="sans-serif" fill="#172c3c">Neutral-plate wetting susceptibility: two independent methods</text>')
    panels.append('<text x="66" y="55" font-size="12" font-family="sans-serif" fill="#536573">Exact Jacobi series versus volume-constrained P1 finite elements; model at assumed sidewall angle 32°</text>')
    panels.append(f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#263746" stroke-width="1.5"/>')
    panels.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#263746" stroke-width="1.5"/>')
    def ycoord(error):
        return bottom-110*(math.log10(error)+5)
    for power,label in ((-3,"10⁻³"),(-4,"10⁻⁴"),(-5,"10⁻⁵")):
        yy=ycoord(10**power)
        panels.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{right}" y2="{yy:.2f}" stroke="#d8e3e7" stroke-dasharray="4 6"/>')
        panels.append(f'<text x="{left-15}" y="{yy+4:.2f}" text-anchor="end" font-size="13" font-family="sans-serif" fill="#36546a">{label}</text>')
    xlocs=[130,290,450,610]
    for x,intervals in zip(xlocs,(16,32,64,128)):
        panels.append(f'<line x1="{x}" y1="{bottom}" x2="{x}" y2="{bottom+5}" stroke="#263746"/>')
        panels.append(f'<text x="{x}" y="{bottom+24}" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#36546a">{intervals}</text>')
    colors=("#bd5634","#247a96")
    for series,(row,color) in enumerate(zip(data["phases"],colors)):
        points=[(x,ycoord(r["absolute_error_mm2"])) for x,r in zip(xlocs,row["fem"])]
        coords=" ".join(f"{x:.1f},{y:.2f}" for x,y in points)
        panels.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for x,y in points:
            panels.append(f'<circle cx="{x}" cy="{y:.2f}" r="4.2" fill="{color}" stroke="#fff" stroke-width="1"/>')
        phase="−" if row["phase"]==-1 else "+"
        q=row["continuum"]["Q_mm2"]
        legend_x=130+series*300
        panels.append(f'<line x1="{legend_x}" y1="410" x2="{legend_x+26}" y2="410" stroke="{color}" stroke-width="3"/>')
        panels.append(f'<text x="{legend_x+36}" y="414" font-size="13" font-family="sans-serif" fill="#253e50">{phase} phase: Q = {q:.10f} mm²</text>')
    panels.append('<text x="393" y="452" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#263746">Number of planar P1 elements (z elements = half as many)</text>')
    panels.append('<text x="22" y="224" transform="rotate(-90 22 224)" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#263746">Absolute energy-coefficient error (mm²)</text>')
    panels.append('<text x="485" y="91" font-size="13" font-family="sans-serif" fill="#5b7180">≈ 4× error reduction per refinement</text>')
    panels.append('</svg>')
    return "\n".join(panels)+"\n"


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--stdout",action="store_true")
    parser.add_argument("--output",default="figures/plate_susceptibility_convergence.svg")
    args=parser.parse_args()
    if args.stdout:
        print(make_svg(),end="")
    else:
        target=Path(args.output)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(make_svg(),encoding="utf-8")
        print(target)
