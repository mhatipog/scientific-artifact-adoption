"""Regenerate small vector comparison from derived, source-attributed saturation CSVs.

No raw videos needed; no time shift or model parameter is fitted.
"""
from __future__ import annotations
import csv
from pathlib import Path
from html import escape

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent / "empirical" / "waveform_registration" / "reference"


def read_rows(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def path_line(pairs, xfun, yfun):
    return " ".join(("M" if i == 0 else "L") +
                    f"{xfun(t):.2f} {yfun(s):.2f}"
                    for i, (t, s) in enumerate(pairs))


def build_svg(mask_rows, papers):
    width, height = 1240, 580
    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
        '<rect width="1240" height="580" fill="#fbfcfe"/>',
        '<text x="58" y="48" font-family="sans-serif" font-size="25" fill="#173245" font-weight="700">Original TIFF masks vs published pore saturation</text>',
        '<text x="58" y="74" font-family="sans-serif" font-size="14" fill="#506778">No time shift · Constant 120 µm depth is a model assumption</text>',
    ]
    for i, deg in enumerate((60, 120)):
        left, top, w, h = (65 if i == 0 else 675), 135, 490, 320
        X = lambda t: left + w * t / 12.8
        Y = lambda s: top + h * (1 - s)
        ours = sorted([
            (float(r["clip_t_s"]), float(r["mask_lo"]),
             float(r["mask_reference"]), float(r["mask_hi"]))
            for r in mask_rows if int(r["angle_deg"]) == deg
            and 0 <= float(r["clip_t_s"]) <= 12.8
        ])
        pub = sorted((float(r["paper_time_s"]),
                      float(r["published_saturation_fraction"]))
                     for r in papers[deg] if 0 <= float(r["paper_time_s"]) <= 12.8)
        if len(ours) < 100 or len(pub) < 100:
            raise ValueError(f"Missing complete source curves: {deg}")
        pieces.append(f'<text x="{left}" y="107" font-family="sans-serif" font-size="21" fill="#203e50" font-weight="700">{deg}° wettability</text>')
        for pct in (0, .2, .4, .6, .8, 1.):
            yy = Y(pct)
            pieces.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{left+w}" y2="{yy:.2f}" stroke="#e2eaf0"/>')
            pieces.append(f'<text x="{left-13}" y="{yy+4:.2f}" text-anchor="end" font-size="12" font-family="sans-serif" fill="#567">{pct:.0%}</text>')
        for sec in range(0, 13, 2):
            xx = X(sec)
            pieces.append(f'<line x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{top+h}" stroke="#edf2f4"/>')
            pieces.append(f'<text x="{xx:.2f}" y="{top+h+21}" text-anchor="middle" font-size="12" font-family="sans-serif" fill="#567">{sec}</text>')
        band = ([(x[0], x[1]) for x in ours] +
                [(x[0], x[3]) for x in reversed(ours)])
        polygon = " ".join(f"{X(t):.2f},{Y(s):.2f}" for t,s in band)
        pieces.append(f'<polygon points="{polygon}" fill="#27a99a" opacity=".18"/>')
        pieces.append(f'<path d="{escape(path_line(pub,X,Y))}" fill="none" stroke="#35495a" stroke-width="2.3"/>')
        pieces.append(f'<path d="{escape(path_line([(x[0],x[2]) for x in ours],X,Y))}" fill="none" stroke="#008d86" stroke-width="2.6"/>')
        pieces.append(f'<text x="{left+w/2}" y="{top+h+48}" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#263f50">Time in image clip (s)</text>')
    pieces += [
        '<line x1="74" y1="548" x2="110" y2="548" stroke="#35495a" stroke-width="3"/>',
        '<text x="117" y="553" font-family="sans-serif" font-size="13">Yan et al. Fig. 4, digitized</text>',
        '<line x1="335" y1="548" x2="371" y2="548" stroke="#008d86" stroke-width="3"/>',
        '<text x="378" y="553" font-family="sans-serif" font-size="13">Independent TIFF mask</text>',
        '<rect x="650" y="540" width="29" height="14" fill="#27a99a" opacity=".25"/>',
        '<text x="687" y="553" font-family="sans-serif" font-size="13">Threshold sensitivity</text>',
        '<text x="1230" y="572" text-anchor="end" font-family="sans-serif" font-size="10" fill="#607681">Yan et al. 2025 CC BY 4.0 · DaRUS-4683</text>',
        '</svg>'
    ]
    return "\n".join(pieces)+"\n"


def main():
    rows = read_rows(HERE / "derived" / "segmented_saturation_10hz.csv")
    papers = {deg:read_rows(REFERENCE / f"fig4_{deg}deg_digitized.csv")
              for deg in (60, 120)}
    output = HERE / "figures" / "phase_saturation_vs_published.svg"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_svg(rows,papers))
    print(f"Regenerated {output}, no fitted timing.")


if __name__ == "__main__":
    main()
