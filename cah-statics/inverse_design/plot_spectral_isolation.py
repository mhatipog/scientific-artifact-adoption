"""Reproduce the SVG showing quartic reversal as the Jacobi spectral gap closes.

Run from this directory: python plot_spectral_isolation.py
Outputs SVG to stdout; no files are modified.
"""
from __future__ import annotations
import math
from spectral_isolation import Channel, record


ORDERS = (15, 25, 51, 75, 101, 151, 201)


def svg_string() -> str:
    rows = [record(Channel(high_order=m,theta_deg=37.1451464361036,
                           high_scale_mm=.046875)) for m in ORDERS]
    x = lambda m: 145 + 800*math.log(m/15)/math.log(201/15)
    top_y = lambda v: 220-.115*v
    bottom_y = lambda value: 410+42.5*(math.log10(.01)-math.log10(value))
    s = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="650" viewBox="0 0 1100 650">',
        '<rect width="1100" height="650" fill="#ffffff"/>',
        '<text x="66" y="41" font-size="25" font-family="sans-serif" font-weight="bold" fill="#182438">Same linear spectrum limit, opposite nonlinear pitchfork classes</text>',
        '<text x="67" y="65" font-size="14" font-family="sans-serif" fill="#4a5b70">Neutral-plate Young-Laplace cylinders: signed high harmonics with width amplitude 0.046875/M^4 mm</text>',
        '<rect x="120" y="89" width="853" height="260" fill="#f8fafd" stroke="#d8e0e8"/>',
        '<rect x="120" y="401" width="853" height="180" fill="#f8fafd" stroke="#d8e0e8"/>',
    ]
    for val in (-1000,-500,0,500,1000):
        y=top_y(val)
        s.append(f'<line x1="120" y1="{y:.2f}" x2="973" y2="{y:.2f}" stroke="#dfe5ec" stroke-width="1"/>')
        s.append(f'<text x="106" y="{y+5:.2f}" text-anchor="end" font-size="13" fill="#4a5b70">{val}</text>')
    for expo in (-2,-3,-4,-5,-6):
        val=10.**expo; y=bottom_y(val)
        s.append(f'<line x1="120" y1="{y:.2f}" x2="973" y2="{y:.2f}" stroke="#dfe5ec"/>')
        s.append(f'<text x="106" y="{y+4:.2f}" text-anchor="end" font-size="13" fill="#4a5b70">1e{expo}</text>')
    for m in (15,51,201):
        xx=x(m)
        s.append(f'<line x1="{xx:.2f}" y1="89" x2="{xx:.2f}" y2="581" stroke="#e5e9ef" stroke-dasharray="4 5"/>')
        s.append(f'<text x="{xx:.2f}" y="599" text-anchor="middle" font-size="14" fill="#25374e">{m}</text>')
    series=[("minus","#22589b"),("plus","#d37521")]
    for phase,color in series:
        points=" ".join(f'{x(row["m"]):.2f},{top_y(row[phase]["B4"]):.2f}' for row in rows)
        s.append(f'<polyline points="{points}" stroke="{color}" stroke-width="2.6" fill="none" stroke-linecap="round" stroke-linejoin="round"/>')
        for row in rows:
            xx=x(row["m"]); y=top_y(row[phase]["B4"])
            s.append(f'<circle cx="{xx:.2f}" cy="{y:.2f}" r="4.6" fill="{color}" stroke="#fff" stroke-width="1"/>')
    points=" ".join(f'{x(row["m"]):.2f},{bottom_y(row["lateral_gap_difference_per_mm2"]):.2f}' for row in rows)
    s.append(f'<polyline points="{points}" stroke="#158779" stroke-width="2.6" fill="none" stroke-linejoin="round"/>')
    for row in rows:
        xx=x(row["m"]);yy=bottom_y(row["lateral_gap_difference_per_mm2"])
        s.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="4.5" fill="#158779" stroke="#fff" stroke-width="1"/>')
    s.extend([
        '<text x="124" y="107" font-size="14" font-family="sans-serif" font-weight="bold" fill="#25374e">Fixed-volume quartic B4 (reduced model)</text>',
        '<line x1="697" y1="106" x2="722" y2="106" stroke="#22589b" stroke-width="3"/>',
        '<text x="727" y="111" font-size="13" font-family="sans-serif" fill="#22589b">minus phase</text>',
        '<line x1="833" y1="106" x2="858" y2="106" stroke="#d37521" stroke-width="3"/>',
        '<text x="863" y="111" font-size="13" font-family="sans-serif" fill="#d37521">plus phase</text>',
        '<text x="124" y="420" font-size="14" font-family="sans-serif" font-weight="bold" fill="#25374e">Difference of first nonzero Jacobi eigenvalues (mm^-2; log scale)</text>',
        '<text x="547" y="620" text-anchor="middle" font-size="15" font-family="sans-serif" fill="#25374e">Odd high-harmonic order M (logarithmic horizontal scale)</text>',
        '<text x="67" y="641" font-size="12" font-family="sans-serif" fill="#6b7788">Mathematical continuum limit only; sub-nanometre high-order widths are not manufacturable.</text>',
        '</svg>',
    ])
    return "\n".join(s) + "\n"


if __name__ == "__main__":
    print(svg_string(),end="")
