"""Reproducible rounded-entry capillary baseline + raw pressure audit (Yan 2025).

Requires no third-party libraries. The wall fillet is an IDEAL 20 um
quarter-circle, not the measured wall polygon. P=water minus air.
This is a 2D principal-curvature construction with a rectangular-depth
curvature closure; NOT a fully resolved 3D free-surface solver.
"""
import csv,json,math
from pathlib import Path
OUT=Path(__file__).parent/"outputs"/"empirical"
OUT.mkdir(parents=True,exist_ok=True)
GAMMA=0.072                 # approximate water-air surface tension N/m
THROAT=150e-6              # designed mask geometry m
PORE=800e-6                # designed pore-body width m
DEPTH=120e-6               # reported microchannel depth m
RADIUS=20e-6               # SEM approximately, not independent exact trace m

def profile(side_deg,plate_deg=None,radius=RADIUS,n=2001):
    if plate_deg is None: plate_deg=side_deg
    side=math.radians(side_deg)
    plate=math.radians(plate_deg)
    rows=[]
    for i in range(n):
        phi=math.pi/2*i/(n-1)
        x=radius*math.sin(phi)
        w=THROAT+2*radius*(1-math.cos(phi))
        p=-GAMMA*(2*math.cos(side+phi)/w+
                  2*math.cos(plate)/DEPTH)
        rows.append(dict(phi_deg=math.degrees(phi),
                    axial_um=x*1e6,
                    nominal_width_um=w*1e6,
                    capillary_kpa=p/1000))
    return rows

def baseline(angle):
    rows=profile(angle)
    peak=max(rows,key=lambda x:x["capillary_kpa"])
    entry=rows[0]
    body=-2*GAMMA*math.cos(math.radians(angle))*(1/PORE+1/DEPTH)/1000
    crossing=(90-angle) if angle<90 else None
    return dict(static_contact_angle_deg=angle,
        throat_capillary_kpa=entry["capillary_kpa"],
        fillet_peak_kpa=peak["capillary_kpa"],
        fillet_peak_axial_um=peak["axial_um"],
        rise_to_peak_kpa=peak["capillary_kpa"]-entry["capillary_kpa"],
        hypothetical_uniform_pore_pressure_kpa=body,
        planar_curvature_reversal_x_um=(RADIUS*math.sin(math.radians(crossing))*1e6 if crossing is not None else None),
        status="MODEL ONLY; advancing wall angle and 3D meniscus are not measured")

def candidate_events(data):
    out=[]
    for f in data["files"]:
        deg=int(f["name"].split("_")[0])
        trace=f.get("decimated_trace",[])
        if deg==60:
            out.append(dict(static_angle_deg=deg,
                 candidate_status="ambiguous: several large transients, exclude from event scoring"))
            continue
        assert len(trace)>100
        i=min(range(len(trace)-1),key=lambda k:trace[k+1][1]-trace[k][1])
        pre=trace[max(0,i-3):i+1]
        post=trace[i+1:min(len(trace),i+5)]
        before=max(pre,key=lambda z:z[1])
        after=min(post,key=lambda z:z[1])
        out.append(dict(static_angle_deg=deg,
            candidate_status="pressure event only, NOT video aligned",
            leading_edge_time_s=trace[i+1][0],
            before_time_s=before[0],after_time_s=after[0],
            before_kpa=before[1],after_kpa=after[1],
            signal_drop_kpa=before[1]-after[1]))
    return out

def main():
    b=[baseline(a) for a in [60,95,120]]
    assert abs(b[0]["planar_curvature_reversal_x_um"]-10)<1e-10
    assert b[1]["planar_curvature_reversal_x_um"] is None
    assert b[2]["rise_to_peak_kpa"]>0 and b[1]["rise_to_peak_kpa"]>0
    with (OUT/"yan_rounded_entry_model_only.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(b[0]))
        writer.writeheader();writer.writerows(b)
    p=OUT/"darus_lowflow_pressure_audit.json"
    observed=candidate_events(json.loads(p.read_text())) if p.exists() else []
    result=dict(method="ideal quarter-circle rounded entrance; no experimental fit",
                source="Yan et al JCIS 2025; DaRUS 10.18419/DARUS-4683",
                nominal_radius_um=20,roughness_um=5,
                instrument_reproducibility_pa=200,
                angle_caveat="paper reports distinct static and advancing angles; static was used here ONLY as sensitivity hypothesis",
                pressure_convention="water minus air, capillary-only; raw pressure has sensor offset and other losses",
                baseline=b,
                observed_candidate_changes=observed,
                comparison_status="NOT validated; model peak-to-throat pressure rise and sensor pre/post change describe different interface locations. Do NOT compute MAE.")
    (OUT/"yan_geometry_vs_pressure_no_fit.json").write_text(json.dumps(result,indent=2))
    # Compact explanatory SVG, no external plot libraries.
    palette=["#158496","#d87345","#7760a6"]
    chunks=['<svg xmlns="http://www.w3.org/2000/svg" width="1250" height="520" viewBox="0 0 1250 520">',
        '<rect width="1250" height="520" fill="#fbfdff"/>',
        '<text x="60" y="42" font-size="25" font-family="Arial" font-weight="700" fill="#172d44">Rounded pore entry | geometric pressure barriers versus sensor changes</text>',
        '<text x="60" y="69" font-size="14" font-family="Arial" fill="#5d6c7b">Calculated barriers use static-angle HYPOTHESES, nominal 20 µm fillet and a depth-curvature closure. Not a parameter-free 3D prediction.</text>']
    def proj(x,y):
        return 100+(x/20)*660, 355-(y+1.4)/3.2*245
    for y in [-1,0,1]:
        Y=proj(0,y)[1]
        chunks.append(f'<line x1="100" x2="760" y1="{Y:.1f}" y2="{Y:.1f}" stroke="#e1e8ed"/>')
        chunks.append(f'<text x="93" y="{Y+4:.1f}" text-anchor="end" font-family="Arial" font-size="12">{y} kPa</text>')
    for ang,col in zip((60,95,120),palette):
        xy=profile(ang,n=151)
        d=" ".join(("M" if i==0 else "L")+f'{proj(v["axial_um"],v["capillary_kpa"])[0]:.2f},{proj(v["axial_um"],v["capillary_kpa"])[1]:.2f}' for i,v in enumerate(xy))
        chunks.append(f'<path d="{d}" stroke="{col}" stroke-width="3" fill="none"/>')
        chunks.append(f'<text x="{770}" y="{145+(ang//35)*41}" fill="{col}" font-size="16" font-family="Arial">{ang}°</text>')
    for x in [0,5,10,15,20]:
        chunks.append(f'<text x="{proj(x,0)[0]:.1f}" y="380" text-anchor="middle" font-family="Arial" font-size="12">{x}</text>')
    chunks.append('<text x="320" y="407" font-size="14" font-family="Arial">Assumed fillet axial coordinate (µm)</text>')
    chunks.append('<text x="843" y="148" font-weight="700" font-size="16" font-family="Arial">Actual pressure signals</text>')
    for i,o in enumerate(observed):
        ang=o["static_angle_deg"]
        if ang==60:
            label="60°  unaligned; no event assigned"
        else:
            label=f'{ang}°  {o["signal_drop_kpa"]:.3f} kPa candidate decrease'
        chunks.append(f'<text x="845" y="{195+i*46}" font-size="15" font-family="Arial" fill="#2c3a4c">{label}</text>')
    chunks.append('<text x="80" y="478" font-size="15" font-family="Arial" font-weight="700" fill="#a44635">Do not compare model rise and measured drop as if they were the same observable.</text>')
    chunks.append('<text x="80" y="499" font-size="12" font-family="Arial" fill="#657488">Experiment: Yan et al. (2025), DOI 10.1016/j.jcis.2025.137884; data: CC BY 4.0, DOI 10.18419/DARUS-4683.</text>')
    chunks.append('</svg>')
    (OUT/"yan_entry_baseline_vs_raw_observables.svg").write_text("\n".join(chunks))
    print("PASS: nominal single-front model, no fitted pressure coefficients; static/advancing caveat")
    print("MODEL",json.dumps(b))
    print("MEASURED CANDIDATES",json.dumps(observed))
    print("NO VALIDATION claimed: location alignment unavailable.")

if __name__=="__main__":
    main()
