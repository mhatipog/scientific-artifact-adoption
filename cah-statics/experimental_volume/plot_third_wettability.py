"""Replot sparse 95-degree optical phase and center-front traces as SVG.

Only standard-library modules are needed. This is a figure of image-derived
metrics; neither plotted curve is a calibrated 3D pore-volume fraction.
"""
import csv
from pathlib import Path

BASE=Path(__file__).parent
DATA=BASE/"derived/third_wettability_95_even_frame_audit.csv"
OUTPUT=BASE/"figures/third_wettability_95_burst.svg"
FRAME_DT=.033328503370285034


def make_svg(rows):
    vals=[dict(row) for row in rows]
    if len(vals)!=24 or [int(q["frame"]) for q in vals]!=list(range(4800,4848,2)):
        raise ValueError("Expected exactly 24 frozen even-frame samples")
    for q in vals:
        f=int(q["frame"])
        if abs(float(q["clip_t_s"])-f*FRAME_DT)>1.1e-6:
            raise ValueError("Unexpected original frame/clip time")
    W,H=900,480
    left,right,top,bottom=96,840,58,380
    def X(t):return left+(t-159.9)/(161.53-159.9)*(right-left)
    def Y(v):return bottom-v*(bottom-top)
    def points(column,normalize=None):
        coords=[]
        for q in vals:
            v=q[column]
            if v is None or v=="":continue
            v=float(v)
            if normalize:v=(v-116)/322
            coords.append(f"{X(float(q['clip_t_s'])):.2f},{Y(v):.2f}")
        return " ".join(coords)
    ticks=[]
    for t in [160.0,160.3,160.6,160.9,161.2,161.5]:
        xx=X(t)
        ticks.append(f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{bottom}" stroke="#e5e7eb"/>'
                     f'<text x="{xx:.1f}" y="404" text-anchor="middle" font-size="14">{t:.1f}</text>')
    for y in [0,.2,.4,.6,.8,1]:
        yy=Y(y)
        ticks.append(f'<line x1="{left}" y1="{yy:.1f}" x2="{right}" y2="{yy:.1f}" stroke="#e5e7eb"/>'
                     f'<text x="84" y="{yy+5:.1f}" text-anchor="end" font-size="14">{y:.1f}</text>')
    t10=160.5966786630411;t90=160.78670406320302
    shade=(f'<rect x="{X(t10):.2f}" y="{top}" width="{X(t90)-X(t10):.2f}" '
           f'height="{bottom-top}" fill="#fff1e6" opacity="0.75"/>')
    out=f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<rect width="{W}" height="{H}" fill="white"/>
<text x="96" y="27" font-family="sans-serif" font-size="21" font-weight="bold" fill="#111827">Single-front pore filling: original 95° candidate TIFF</text>
<g font-family="sans-serif" fill="#374151">
{shade}
{''.join(ticks)}
<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#374151" stroke-width="1.5"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#374151" stroke-width="1.5"/>
<polyline fill="none" stroke="#2764a5" stroke-width="3.2" stroke-linejoin="round" points="{points('mask_reference')}"/>
<polyline fill="none" stroke="#ce6330" stroke-width="2.7" stroke-dasharray="7 4" stroke-linejoin="round" points="{points('front_x50_px',True)}"/>
<text x="468" y="438" text-anchor="middle" font-size="16">Time within archived TIFF clip (s)</text>
<text x="18" y="214" text-anchor="middle" transform="rotate(-90 18 214)" font-size="15">2D optical metric (dimensionless)</text>
<line x1="525" y1="50" x2="554" y2="50" stroke="#2764a5" stroke-width="3"/>
<text x="560" y="54" font-size="12">Phase mask / pore image area</text>
<line x1="711" y1="50" x2="739" y2="50" stroke="#ce6330" stroke-width="2.8" stroke-dasharray="7 4"/>
<text x="746" y="54" font-size="12">Front extent</text>
<text x="98" y="463" font-size="12">Raw DaRUS TIFF 351440 · 24 measured frames · uniform-depth volume NOT established</text>
<text x="611" y="76" font-size="12" fill="#9a3412">10–90% optical change: 0.190 s</text>
</g></svg>'''
    return out


def main():
    with DATA.open(newline="") as f:svg=make_svg(csv.DictReader(f))
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(svg+"\n",encoding="utf8")
    print(OUTPUT)


if __name__=="__main__":
    main()
