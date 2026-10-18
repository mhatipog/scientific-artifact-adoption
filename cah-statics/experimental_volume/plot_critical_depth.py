"""Offline, dependency-free SVG of Round-13 conditional depth thresholds.

Usage:
    python cah-statics/experimental_volume/plot_critical_depth.py
Reads frozen derived CSV; does not download original TIFF or invoke Actions.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape

DATA = Path(__file__).parent / "derived" / "critical_depth_95_windows.csv"
OUT = Path(__file__).parent / "figures" / "critical_depth_95_guarded.svg"


def svg_string(rows):
    rows = list(rows)
    tags = [(100, "0.4"), (180, "0.4"), (180, "0.5"),
            (180, "0.6"), (260, "0.6")]
    lookup = {(int(r["frame_start"]), int(r["count_threshold"]),
               str(r["fraction_threshold"])):r for r in rows}
    if len(rows)!=10 or len(lookup)!=10:
        raise ValueError("Expected exactly ten distinct measured mask/window cases")
    w,h=960,520
    left,right,top,bottom=90,910,130,393
    lo,hi=.5,20.
    def y(value):
        return bottom-(math.log10(value)-math.log10(lo)) / (
            math.log10(hi)-math.log10(lo))*(bottom-top)
    def t(x, yy, label, size=14, color="#273247", anchor="middle",
          weight="normal"):
        return (f'<text x="{x}" y="{yy}" fill="{color}" font-size="{size}" '
                f'font-family="DejaVu Sans,Arial,sans-serif" '
                f'font-weight="{weight}" text-anchor="{anchor}">'
                +escape(str(label))+'</text>')
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
         f'viewBox="0 0 {w} {h}" role="img" '
         f'aria-label="Lost-pixel-guarded critical thicknesses for five optical masks">'
         '<rect width="960" height="520" fill="#ffffff"/>']
    svg += [t(38,38,"95° pore invasion: pump-only depth certificates",23,
              "#152139","start","bold"),
            t(38,67,"Original frames 4820–4822 versus 4822–4824; each 66.7 ms",15,
              "#586578","start")]
    for v in (.5,1,2,5,10,20):
        yy=y(v)
        svg.append(f'<line x1="{left}" x2="{right}" y1="{yy:.2f}" '
                   f'y2="{yy:.2f}" stroke="#d7dce4" stroke-width="1"/>')
        svg.append(t(left-12,round(yy+5,2),str(v),13,"#5b6676","end"))
    svg.append(t(36,105,"Critical added depth (µm; logarithmic scale)",13,
                 "#42516b","start"))
    colors={4820:"#ce8453",4822:"#28888d"}
    for idx,(counts,frac) in enumerate(tags):
        x=163+idx*151
        for start,offset in ((4820,-26),(4822,13)):
            row=lookup[(start,counts,frac)]
            value=float(row["lost_pixel_guarded_critical_um"])
            yy=y(value); xx=x+offset
            svg.append(f'<rect x="{xx}" y="{yy:.2f}" width="31" '
                       f'height="{bottom-yy:.2f}" rx="2" fill="{colors[start]}"/>')
            svg.append(t(xx+15.5,round(yy-7,2),f'{value:.2f}',12,
                         colors[start],weight="bold"))
        svg.append(t(x+4,420,f'{counts} / {frac}',14,"#273247"))
    ref=y(2.)
    svg.append(f'<line x1="{left}" x2="{right}" y1="{ref:.2f}" '
               f'y2="{ref:.2f}" stroke="#56516d" stroke-width="2" '
               'stroke-dasharray="7 5"/>')
    svg += [t(107,115,"2 µm is a HYPOTHETICAL incremental-depth target, not measured",
              11,"#675d7e","start"),
            f'<rect x="90" y="449" width="13" height="13" fill="{colors[4820]}"/>',
            t(112,460,"4820–4822: mask reversals up to 10.6%",13,
              "#45556b","start"),
            f'<rect x="495" y="449" width="13" height="13" fill="{colors[4822]}"/>',
            t(518,460,"4822–4824: reversals below 0.16%",13,
              "#45556b","start"),
            t(91,493,"Guard credits disappearing pixels with releasing ALL 120 µm; retained-pixel thinning remains unobserved.",
              12,"#677389","start"),
            '</svg>']
    return '\n'.join(svg)+'\n'


def main():
    with DATA.open(newline="") as stream:
        data=list(csv.DictReader(stream))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(svg_string(data),encoding="utf-8")
    print(str(OUT))


if __name__=="__main__":
    main()
