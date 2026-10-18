"""Offline SVG of the nominal 0.17-µm fifth-order isocorrelated pair."""
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
from phase_pair_window import Geometry, pair

OUT=Path(__file__).parent/"figures"/"phase_pair_robust_window.svg"


def svg_string():
    w,h=1020,560
    blue="#188895";orange="#cb7951";dark="#22304b";muted="#69778b"
    shapes=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img" '
            'aria-label="Opposite local quartic coefficients for nearly identical roughness channels">'
            f'<rect x="0" y="0" width="{w}" height="{h}" fill="#ffffff"/>']
    def txt(x,y,s,size=14,color=dark,anchor="middle",weight="normal"):
        return f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" ' +(
            f'font-family="DejaVu Sans,Arial,sans-serif" font-size="{size}" '
            f'font-weight="{weight}" fill="{color}">{escape(str(s))}</text>')
    shapes.append(txt(36,39,"Fifth-order-identical roughness, opposite local onset classes",
                      23,dark,"start","bold"))
    shapes.append(txt(36,67,"Planar Young–Laplace design; θY = 32° nominal and UNMEASURED",
                      15,muted,"start"))
    # Panel 1: deviation in total channel width from its shared fundamental.
    shapes.append(txt(270,112,"High-harmonic width correction (µm)",16,dark,weight="bold"))
    u=np.linspace(0,2*np.pi,801)
    x=80+400*u/(2*np.pi)
    ripple=.17*(np.cos(5*u)-np.cos(7*u))
    yscale=lambda v:340-v*310
    for tick in (-.3,0,.3):
        yy=yscale(tick)
        shapes.append(f'<line x1="80" y1="{yy:.2f}" x2="480" y2="{yy:.2f}" '
                      'stroke="#e5e9ef" stroke-width="1"/>')
        shapes.append(txt(70,yy+5,str(tick),12,muted,"end"))
    for sign,color in ((-1,orange),(1,blue)):
        yy=yscale(sign*ripple)
        p=" ".join(f'{a:.2f},{b:.2f}' for a,b in zip(x,yy))
        shapes.append(f'<polyline fill="none" stroke="{color}" '
                      f'stroke-width="2.2" points="{p}"/>')
    for ix,label in ((80,"0"),(280,"150"),(480,"300")):
        shapes.append(txt(ix,370,label,13,muted))
    shapes.append(txt(280,393,"Axial distance in one 300 µm period",13,muted))
    shapes.append(txt(280,425,"Max paired width difference: 0.657 µm",13,dark))
    # Panel 2: independent analytic B4 at five intrinsic sidewall angles.
    shapes.append(txt(760,112,"Quartic coefficient B₄ (model units)",16,dark,weight="bold"))
    vals=[]
    for angle in (30,31,32,33,34):
        r=pair(Geometry(theta_deg=angle),.00017,.00017)
        if not r["ok"]:
            raise ValueError(f"Nominal pair has no valid local onset at {angle}")
        vals.append((angle,r["B4_minus"],r["B4_plus"]))
    axis_x=lambda theta:600+(theta-30)*80
    axis_y=lambda val:240-val*110/4000
    for tick in (-4000,-2000,0,2000,4000):
        yy=axis_y(tick)
        shapes.append(f'<line x1="600" y1="{yy:.2f}" x2="920" y2="{yy:.2f}" '
                      f'stroke="{"#9ba7b8" if tick==0 else "#e5e9ef"}" '
                      'stroke-width="1"/>')
        shapes.append(txt(590,yy+4,str(tick),12,muted,"end"))
    for column,color in ((1,orange),(2,blue)):
        points=" ".join(f'{axis_x(v[0]):.2f},{axis_y(v[column]):.2f}'
                        for v in vals)
        shapes.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.4" '
                      f'points="{points}"/>')
        for v in vals:
            shapes.append(f'<circle cx="{axis_x(v[0]):.2f}" '
                          f'cy="{axis_y(v[column]):.2f}" r="3.8" fill="{color}"/>')
    for angle in (30,31,32,33,34):
        shapes.append(txt(axis_x(angle),370,str(angle),13,muted))
    shapes.append(txt(760,393,"Assumed sidewall Young angle (degrees)",13,muted))
    shapes.append(txt(760,425,"Local, not a finite-jump or 3D wetting prediction",12,muted))
    # Shared legend and explicit no-experiment claim.
    shapes.append(f'<line x1="145" y1="458" x2="175" y2="458" stroke="{orange}" stroke-width="3"/>')
    shapes.append(txt(184,463,"Phase −",14,dark,"start"))
    shapes.append(f'<line x1="297" y1="458" x2="327" y2="458" stroke="{blue}" stroke-width="3"/>')
    shapes.append(txt(335,463,"Phase +",14,dark,"start"))
    shapes.append(txt(36,500,"Exact: same PSD, extrema and correlations through order five. Sixth order differs.",
                      13,muted,"start"))
    shapes.append(txt(36,525,"Numerical: 405 sampled uncertainty cases pass; nearest fold-creation headroom ~38 nm.",
                      13,muted,"start"))
    shapes.append('</svg>')
    return '\n'.join(shapes)+'\n'


def main():
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(svg_string(),encoding="utf-8")
    print(OUT)


if __name__=="__main__":
    main()
