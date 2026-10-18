"""Render the Round-18 local nonlinear KKT convergence data as standalone SVG."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

BASE=Path(__file__).resolve().parent


def make_svg(rows: list[dict]) -> str:
    grouped={-1:[],1:[]}
    for row in rows:
        grouped[int(row["phase"])].append(row)
    for phase in grouped:
        grouped[phase].sort(key=lambda r:int(r["nt"]))
        if len(grouped[phase])!=5:
            raise ValueError("Need the five declared diagonal meshes for each phase")
    xloc=[130,270,410,550,690]
    def ycoord(v: float):
        return 365-160*(math.log10(v)+.7)
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="850" height="480" viewBox="0 0 850 480">',
    '<rect width="850" height="480" fill="#ffffff"/>',
    '<text x="65" y="34" font-size="20" font-family="sans-serif" fill="#193242">Nonlinear 3D Young–Laplace: independent benchmark convergence</text>',
    '<text x="65" y="57" font-size="13" font-family="sans-serif" fill="#526775">Error of finite-angle, fixed-volume Q against the exact neutral-plate continuum Jacobi result</text>',
    '<line x1="104" y1="85" x2="104" y2="365" stroke="#263746" stroke-width="1.5"/>',
    '<line x1="104" y1="365" x2="730" y2="365" stroke="#263746" stroke-width="1.5"/>']
    for p in (.5,1,2,5,10):
        y=ycoord(p)
        svg.append(f'<line x1="104" y1="{y:.2f}" x2="730" y2="{y:.2f}" stroke="#dbe5e8" stroke-dasharray="3 5"/>')
        svg.append(f'<text x="92" y="{y+5:.2f}" text-anchor="end" font-size="13" font-family="sans-serif" fill="#36546a">{p:g}%</text>')
    names={-1:"Minus Fourier phase",1:"Plus Fourier phase"}
    colors={-1:"#bb553d",1:"#227d95"}
    for phase in (-1,1):
        rr=grouped[phase]
        points=[(x,ycoord(float(v["relative_Q_error_pct"]))) for x,v in zip(xloc,rr)]
        svg.append('<polyline points="'+ " ".join(f"{x},{y:.2f}" for x,y in points)+f'" fill="none" stroke="{colors[phase]}" stroke-width="2.7"/>')
        for x,y in points:
            svg.append(f'<circle cx="{x}" cy="{y:.2f}" r="4.7" fill="{colors[phase]}" stroke="white" stroke-width="1"/>')
        lx=170 if phase==-1 else 495
        svg.append(f'<line x1="{lx}" y1="399" x2="{lx+24}" y2="399" stroke="{colors[phase]}" stroke-width="3"/>')
        svg.append(f'<text x="{lx+33}" y="403" font-size="13" font-family="sans-serif" fill="#36546a">{names[phase]}</text>')
    for x,lab in zip(xloc,("65×5","97×7","129×9","193×13","257×17")):
        svg.append(f'<text x="{x}" y="386" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#36546a">{lab}</text>')
    svg.extend([
     '<text x="412" y="442" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#263746">Planar × depth nodes; both directions refined together</text>',
     '<text x="21" y="220" transform="rotate(-90 21 220)" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#263746">Absolute relative Q error (logarithmic scale)</text>',
     '<text x="412" y="461" text-anchor="middle" font-size="12" font-family="sans-serif" fill="#687d88">89.9° finite-angle Q; analytical reference at 90°. Not a nonlinear bifurcation-class test.</text>',
     '</svg>'])
    return "\n".join(svg)+"\n"


def read_rows() -> list[dict]:
    with (BASE/"derived/nonlinear3d_local_diagonal.csv").open(newline="",encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--stdout",action="store_true")
    parser.add_argument("--output",default=str(BASE/"figures/nonlinear3d_local_convergence.svg"))
    args=parser.parse_args()
    svg=make_svg(read_rows())
    if args.stdout:
        print(svg,end="")
    else:
        target=Path(args.output)
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(svg,encoding="utf-8")
        print(target)
