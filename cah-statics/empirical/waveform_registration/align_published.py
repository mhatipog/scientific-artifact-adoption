"""Digitize Yan et al. 2025 Figure 5 blue curves and match ORIGINAL FISO traces.

No theoretical model or jump parameter is used in the registration.
Only the time offset is optimized: no gain or vertical pressure offset.
All scales are approximate paper-figure axes; retain audit residuals.
Source: Yan et al. JCIS 696 (2025), Fig. 5, CC BY 4.0;
raw pressure: DaRUS 10.18419/DARUS-4683, CC BY 4.0.
"""
from pathlib import Path
from datetime import datetime
from urllib.request import urlopen,Request
import hashlib,json,io,csv
import numpy as np
from PIL import Image
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE/"outputs"
OUT.mkdir(exist_ok=True)
IMG=OUT/"published_fig5_source_page.png"
FILEIDS={60:351442,95:351431,120:351436}
# Explicit calibrated axes (pixels at figure_source.py 1.8x PDF rendering).
# Box x0,y0,x1,y1; graph axes correspond to (t0,pressure0) and slopes.
AXES={
  60:dict(box=(120,112,329,299),x0=119,px_per_second=18.75,y0=109,p0_pa=-100,dy_px=194,dp_pa=-600),
  95:dict(box=(123,365,338,550),x0=122,px_per_second=14.55,y0=556,p0_pa=400,dy_px=-166,dp_pa=600),
 120:dict(box=(123,620,339,803),x0=122,px_per_second=24.40,y0=808,p0_pa=500,dy_px=-193,dp_pa=700)
}
USERAGENT={"User-Agent":"IndependentScientificReferenceMatcher/1.0"}

def published_trace(deg,pixels):
    c=AXES[deg];x1,y1,x2,y2=c["box"]
    r=pixels[y1:y2,x1:x2,:].astype(np.int16)
    mask=(r[:,:,2]-r[:,:,0]>45)&(r[:,:,2]-r[:,:,1]>30)&(r[:,:,2]>140)
    xs=np.arange(x1,x2)
    yy=np.array([(y1+np.median(np.flatnonzero(mask[:,j]))) if mask[:,j].any()
                 else np.nan for j in range(x2-x1)])
    known=np.isfinite(yy)
    if known.mean()<.96:raise RuntimeError(f"insufficient blue trace pixels: {deg} {known.mean():.3f}")
    yy=np.interp(xs,xs[known],yy[known])
    yy=gaussian_filter1d(median_filter(yy,3),1.8)
    tt=(xs-c["x0"])/c["px_per_second"]
    pp=(c["p0_pa"]+(yy-c["y0"])*c["dp_pa"]/c["dy_px"])/1000
    keep=tt>=0
    return tt[keep],pp[keep]

def raw_pressure(deg):
    url=f"https://darus.uni-stuttgart.de/api/access/datafile/{FILEIDS[deg]}"
    with urlopen(Request(url,headers=USERAGENT),timeout=35) as resp:
        raw=resp.read(8_000_001)
    if len(raw)>8_000_000 or b"[Data]" not in raw:raise RuntimeError("Expected original FISO TXT absent/oversized")
    content=raw.decode("utf-8-sig")
    hdr,measure=content.split("[Data]",1)
    if "\nUnit\tkPa" not in hdr:raise RuntimeError("Pressure exported unit is not kPa")
    tstart=None;ts=[];ps=[];invalid=0
    for line in measure.splitlines():
        row=line.split("\t")
        if len(row)<4:continue
        if row[3].strip()!="Valid":invalid+=1;continue
        try:
            dt=datetime.fromisoformat(row[0]+"T"+row[1])
            pressure=float(row[2])
        except ValueError:continue
        if tstart is None:tstart=dt
        ts.append((dt-tstart).total_seconds());ps.append(pressure)
    tt,pp=np.array(ts),np.array(ps)
    if len(tt)<1000 or np.any(np.diff(tt)<=0):raise RuntimeError("Invalid original pressure timestamps")
    # Mean filter reduces 62.5 Hz sensor jitter; source already figure-averaged.
    dt=float(np.median(np.diff(tt)))
    smoothed=gaussian_filter1d(pp,.30/dt)
    return tt,smoothed,{"original_md5":hashlib.md5(raw).hexdigest(),
           "valid_samples":len(tt),"invalid_samples":invalid,
           "native_unit":"kPa","sampling_interval_s":dt}

def match(tpub,ppub,tlogged,plogged):
    offsets=np.arange(0,tlogged[-1]-tpub[-1],.25)
    loss=np.empty_like(offsets)
    for i,shift in enumerate(offsets):
        read=np.interp(shift+tpub,tlogged,plogged)
        loss[i]=np.mean((read-ppub)**2)
    ix=int(np.argmin(loss))
    # Refine near the coarse winner; still time shift ONLY.
    from scipy.optimize import minimize_scalar
    def error(s):
        read=np.interp(s+tpub,tlogged,plogged)
        return np.mean((read-ppub)**2)
    fit=minimize_scalar(error,bounds=(max(0,offsets[ix]-.3),
                        min(tlogged[-1]-tpub[-1],offsets[ix]+.3)),
                        method="bounded",options={"xatol":1e-4})
    # Compare with local minima separated by 20 seconds.
    minima,_=find_peaks(-loss,distance=80)
    closest=sorted(minima,key=lambda k:loss[k])
    # Preserve potential endpoint minima if they ever beat interior candidates
    alt=[(float(offsets[k]),float(np.sqrt(loss[k])))
         for k in closest if abs(offsets[k]-fit.x)>20][:4]
    best=float(fit.x)
    return best,float(np.sqrt(fit.fun)),alt,np.interp(best+tpub,tlogged,plogged)

def main():
    pixels=np.asarray(Image.open(IMG).convert("RGB"))
    report={"status":"publication waveform matched to raw sensor; CAMERA->PUBLICATION time mapping still needs independent check",
            "figure_source":"Yan et al (2025), Fig. 5, CC BY 4.0",
            "dataset":"doi:10.18419/DARUS-4683","pressure_gain_fit":False,
            "pressure_offset_fit":False,"model_parameter_fit":False,
            "matches":[]}
    fig,axes=plt.subplots(3,2,figsize=(14,9),layout="constrained")
    for row,deg in enumerate([60,95,120]):
        t,p=published_trace(deg,pixels)
        tt,pp,provenance=raw_pressure(deg)
        shift,rmse,alternatives,matchraw=match(t,p,tt,pp)
        assert np.isfinite(rmse) and rmse<.065, (deg,shift,rmse)
        result=dict(angle_deg=deg,matching_time_origin_sensor_s=shift,
                    rms_waveform_disagreement_kPa=rmse,
                    next_separate_candidate=alternatives[0] if alternatives else None,
                    number_published_pixels=len(t),duration_figure_s=float(t[-1]),
                    original_sensor=provenance,
                    absolute_camera_synchronization="NOT ESTABLISHED")
        report["matches"].append(result)
        ax=axes[row,0]
        ax.plot(t,p,lw=1.5,label="published Fig. 5, digitized blue waveform")
        ax.plot(t,matchraw,lw=2,label="actual FISO (only time shift)")
        ax.set(title=f"{deg}° | Time origin = sensor {shift:.2f}s; RMS {1000*rmse:.0f} Pa",
               xlabel="time on paper Fig. 5 (s)",ylabel="pressure (kPa)")
        ax.legend(fontsize=7);ax.grid(alpha=.16)
        ax=axes[row,1]
        ax.plot(tt[::55],pp[::55],lw=1,label="original pressure recording")
        ax.axvspan(shift,shift+t[-1],alpha=.18,label="matching published interval")
        ax.set(xlabel="sensor recording time (s)",ylabel="pressure (kPa)")
        ax.grid(alpha=.16);ax.legend(fontsize=8)
        # Figures derived from authors' CC-BY work, not original manuscript.
        with (OUT/f"figure5_{deg}deg_digitized.csv").open("w",newline="") as f:
            writer=csv.writer(f);writer.writerow(["paper_time_s","paper_pressure_kPa","original_sensor_pressure_kPa"])
            writer.writerows(zip(t,p,matchraw))
        print("MATCH",deg,"start",shift,"rmse_kPa",rmse,"alternatives",alternatives[:2],flush=True)
    fig.suptitle("Original published experiment-to-pressure clock matching\nNo capillary-model or pressure-baseline fitting",
                 fontsize=17,weight="bold")
    fig.savefig(OUT/"independently_registered_experimental_pressure.png",dpi=165)
    fig.savefig(OUT/"independently_registered_experimental_pressure.svg")
    plt.close(fig)
    (OUT/"registration_audit.json").write_text(json.dumps(report,indent=2))
    print("PASS matching source experimental waveforms; no video clock or capillary pressure model claimed",flush=True)

if __name__=="__main__":
    main()
