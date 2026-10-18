"""Experimental image-to-pressure time bridge, with zero capillary-model fitting.

Published Fig. 4 (saturation) vs separately extracted TIFF optical response:
ONLY relative experimental clock shift optimized. The Fig. 5 pressure
recording origin is independently identified by align_published.py.
Source: Yan et al., JCIS 696 (2025), CC BY 4.0,
doi:10.1016/j.jcis.2025.137884 and DaRUS-4683.
"""
from pathlib import Path
import csv,json
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE/"outputs"
OUT.mkdir(parents=True,exist_ok=True)
ORIGINAL=(np.asarray(Image.open(OUT/"published_fig4_source_page.png").convert("RGB"),dtype=np.int16)
          if (OUT/"published_fig4_source_page.png").exists() else None)
REGISTER=json.loads((OUT/"registration_audit.json").read_text()
       if (OUT/"registration_audit.json").exists() else
       (HERE/"reference"/"registration_audit.json").read_text())
CLOCK={r["angle_deg"]:r["matching_time_origin_sensor_s"] for r in REGISTER["matches"]}

def publication_saturation(deg):
    offline=HERE/"reference"/f"fig4_{deg}deg_digitized.csv"
    if ORIGINAL is None:
        a=np.genfromtxt(offline,delimiter=",",skip_header=1)
        return a[:,0],a[:,1]
    x1,y1,x2,y2=114,124,360,323
    c=ORIGINAL[y1:y2,x1:x2]
    if deg==60:
        mask=(c[:,:,2]-c[:,:,0]>70)&(c[:,:,2]-c[:,:,1]>60)&(c[:,:,2]>150)
    else:
        mask=(c[:,:,0]>160)&(c[:,:,1]>130)&(c[:,:,2]<130)
    # Exclude the blue/yellow sample markers of the publication's plot legend.
    mask[(270-y1):,(297-x1):]=False
    x=np.arange(x1,x2)
    y=np.array([y1+np.median(np.flatnonzero(mask[:,i])) if mask[:,i].any()
                else np.nan for i in range(x2-x1)])
    ok=np.isfinite(y)
    if ok.sum()<175:raise RuntimeError("Insufficient original experimental dots")
    return (x[ok]-113)/(249/14), (323-y[ok])/200

def main():
    fig,axs=plt.subplots(2,2,figsize=(13,8),layout="constrained")
    with (HERE.parent/"yan_video"/"observed_video_summary.csv").open() as f:
        previous=list(csv.DictReader(f))
    results=[]
    for i,deg in enumerate((60,120)):
        t,s=publication_saturation(deg)
        data=np.genfromtxt(HERE/f"derived_optical_{deg}.csv",delimiter=",",skip_header=1)
        tv,op=data[:,0],data[:,2]
        if len(tv)!=(381 if deg==60 else 658):raise RuntimeError("Original TIFF-derived data incomplete")
        b,fin=(.05,1.) if deg==60 else (0.,.92)
        core=(t>=1)&(t<=(10 if deg==60 else 9))
        tcore,score=t[core],s[core]
        def pred(dt,base=b,final=fin):
            return base+(final-base)*np.clip(np.interp(tcore-dt,tv,op),0,1)
        fit=minimize_scalar(lambda dt:np.mean((score-pred(dt))**2),
                            bounds=(-.6,.6),method="bounded",
                            options={"xatol":1e-5})
        shift=float(fit.x);rms=float(np.sqrt(fit.fun))
        shifts=[]
        for base in ([.03,.05,.07] if deg==60 else [0.,.02,.04]):
            for final in ([.98,1,1.02] if deg==60 else [.90,.92,.94]):
                q=minimize_scalar(lambda dt:np.mean((score-pred(dt,base,final))**2),
                     bounds=(-.6,.6),method="bounded")
                shifts.append(float(q.x))
        paper_pressure=np.genfromtxt(
            (OUT/f"figure5_{deg}deg_digitized.csv") if
            (OUT/f"figure5_{deg}deg_digitized.csv").exists() else
            (HERE/"reference"/f"fig5_{deg}deg_matched_fiso.csv"),
            delimiter=",",skip_header=1)
        tp,p=paper_pressure[:,0],paper_pressure[:,2] # real sensor pressure, figure-aligned
        dpdt=np.gradient(gaussian_filter1d(p,.20/np.median(np.diff(tp))),tp)
        valid=(tp>=4)&(tp<=(9.8 if deg==60 else 8.2))
        largest_pressure_s=float(tp[np.where(valid)[0][np.argmin(dpdt[valid])]])
        doptdt=np.gradient(gaussian_filter1d(op,3.),tv)
        watch=(tv>=5.5)&(tv<=(10 if deg==60 else 9))
        largest_optical_s=float(tv[np.where(watch)[0][np.argmax(doptdt[watch])]]+shift)
        row=next(x for x in previous if str(deg)+"_Q1_" in x["experiment"])
        t10=float(row["t10_s"]);t90=float(row["t90_s"])
        press_diff=lambda sh:float(np.interp(t90+sh,tp,p)-np.interp(t10+sh,tp,p))
        rang=[press_diff(shift+u) for u in (-.4,-.2,0,.2,.4)]
        result={"angle_deg":deg,"publication_time_minus_video_time_s":shift,
                "optical_publication_saturation_rms_percentage_points":100*rms,
                "sensitivity_of_camera_shift_to_brightness_normalization_s":[min(shifts),max(shifts)],
                "pressure_logger_origin_for_publication_time_s":CLOCK[deg],
                "time_of_fastest_pressure_fall_published_s":largest_pressure_s,
                "time_of_fastest_optical_filling_published_s":largest_optical_s,
                "optical_minus_pressure_peak_time_s":largest_optical_s-largest_pressure_s,
                "pressure_difference_during_optical_10_to_90pct_kPa_if_clocks_match":press_diff(shift),
                "time_shift_plusminus_0p4s_pressure_difference_range_kPa":[min(rang),max(rang)],
                "interpretation":"experimental time-axis bridge ONLY; optical intensity is not calibrated pore volume"}
        results.append(result)
        ax=axs[i,0]
        ax.plot(t,s*100,"o",ms=2.3,alpha=.7,label="published experimental saturation")
        ax.plot(tv+shift,100*(b+(fin-b)*np.clip(op,0,1)),lw=2,
                label="original TIFF optical response")
        ax.set(xlim=(0,13),ylim=(-4,105),xlabel="paper experimental time (s)",
               ylabel="saturation / optical proxy (%)",
               title=f"{deg}°: clock shift +{shift:.2f}s; disagreement {100*rms:.1f} pp")
        ax.legend(fontsize=8);ax.grid(alpha=.15)
        ax=axs[i,1]
        ax.plot(tp,p,lw=1.8,label="original recorded inlet pressure")
        ax.axvline(largest_pressure_s,color="#b35253",ls="--",label="fastest pressure fall")
        ax.axvline(largest_optical_s,color="#438e75",ls=":",label="fastest optical invasion")
        ax.set(xlim=(0,11 if deg==60 else 9),xlabel="paper experimental time (s)",
               ylabel="measured inlet pressure (kPa)",
               title=f"{deg}°: optical minus pressure timing {largest_optical_s-largest_pressure_s:+.2f}s")
        ax.legend(fontsize=8);ax.grid(alpha=.15)
        print("BRIDGE",json.dumps(result),flush=True)
    fig.suptitle("Two independent experimental time references: camera and pressure\n"
                 "No fitted capillary-statics parameter; optical intensity is not liquid volume",
                 fontsize=16,weight="bold")
    fig.savefig(OUT/"experimental_camera_sensor_bridge.png",dpi=165)
    fig.savefig(OUT/"experimental_camera_sensor_bridge.svg")
    plt.close(fig)
    report={"source":"Yan et al JCIS (2025), Figs 4–5; DaRUS DOI 10.18419/DARUS-4683",
            "pressure_amplitude_fit":False,"pressure_offset_fit":False,
            "theoretical_parameter_fit":False,
            "estimated_only":"relative alignment of published saturation and optical intensity",
            "camera_clock_caveat":"paper Q1 curves must share experimental time origin; no absolute camera acquisition timestamp in TIFF",
            "results":results}
    (OUT/"experimental_camera_pressure_bridge.json").write_text(json.dumps(report,indent=2))
    assert all(r["optical_publication_saturation_rms_percentage_points"]<7 for r in results)
    assert all(abs(r["publication_time_minus_video_time_s"])<.35 for r in results)
    print("PASS two independent data-only bridges: original FISO->paper Fig5; original TIFF->paper Fig4")
if __name__=="__main__":
    main()
