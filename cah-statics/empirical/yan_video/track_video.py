"""Quantify the observed optical filling time, not a calibrated 3D volume."""
import csv,json,re
from pathlib import Path
import numpy as np,tifffile
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/"outputs"
Q=.01*1000/3600    # Q1: 0.01 mL/h -> uL/s
VPORE=800*800*120/1e9 # nominal 800x800x120 um pore -> uL

def crossing(t,y,value):
    ii=np.where(y>=value)[0]
    if not len(ii):return None
    j=int(ii[0])
    if j==0:return float(t[0])
    i=j-1
    return float(t[i]+(value-y[i])*(t[j]-t[i])/(y[j]-y[i]+1e-12))

def measure(path,chosen):
    with tifffile.TiffFile(path) as tif:
        pages=tif.pages
        first=pages[0].asarray().astype(np.float32)
        last=pages[-1].asarray().astype(np.float32)
        desc=str(pages[0].tags["ImageDescription"].value)
        match=re.search(r"finterval\s*=\s*([0-9.]+)",desc)
        dt=float(match.group(1)) if match else 1/30
        H,W=first.shape
        yy,xx=np.indices(first.shape)
        # Documented fixed interior-pore ROI. Pixel bounds may need refining
        # by independent wall segmentation for publication-grade volumes.
        roi=(xx>=110)&(xx<405)&(yy>=95)&(yy<390)
        D=first-last
        positive=D[(D>0)&roi]
        typical=float(np.median(positive))
        thresholds=[.12,.22,.35]
        masks=[roi&(D>max(40,t*typical)) for t in thresholds]
        if any(m.sum()<20000 for m in masks):raise RuntimeError("Insufficient dynamic pore interior")
        denominator=[float(np.sum(D[m])) for m in masks]
        y=np.zeros((len(pages),3))
        for n,page in enumerate(pages):
            frame=page.asarray().astype(np.float32)
            delta=first-frame
            for j,mask in enumerate(masks):
                y[n,j]=np.sum(delta[mask])/denominator[j]
    t=np.arange(len(y))*dt
    results=[]
    for j,cutoff in enumerate(thresholds):
        env=np.maximum.accumulate(gaussian_filter1d(y[:,j],1.4,mode="nearest"))
        t10=crossing(t,env,.1);t90=crossing(t,env,.9)
        if t10 is None or t90 is None or t90<=t10:raise RuntimeError("No reliable 10–90% filling interval")
        elapsed=t90-t10
        results.append({"mask_threshold":cutoff,"changed_pixels":int(masks[j].sum()),
                        "t10_s":t10,"t90_s":t90,"observed_10_90_duration_s":elapsed,
                        "pump_delivery_in_interval_uL":Q*elapsed,
                        "80_percent_of_nominal_pore_volume_uL":.8*VPORE,
                        "estimated_optical_vs_pump_flow_factor":.8*VPORE/(Q*elapsed)})
    OUT.mkdir(exist_ok=True)
    with (OUT/"framewise_optical_response.csv").open("w",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(["t_since_video_clip_start_s","12pct_mask","22pct_mask","35pct_mask"])
        writer.writerows(zip(t,y[:,0],y[:,1],y[:,2]))
    fig,ax=plt.subplots(figsize=(10,5),layout="constrained")
    for j,cutoff in enumerate(thresholds):
        ax.plot(t,y[:,j],label=f"optical response, {int(100*cutoff)}% mask")
    middle=results[1]
    ax.axvspan(middle["t10_s"],middle["t90_s"],alpha=.13,label="10–90% optical interval")
    baseline=.1+(t-middle["t10_s"])*Q/VPORE
    ax.plot(t,baseline,"k--",lw=1.5,label="nominal pump-only saturation slope")
    ax.set(xlabel="Time in original microscope clip (s)",
           ylabel="Normalized optical pore-filling signal",
           title="Measured pore filling from 381 original microscope frames",
           ylim=(-.10,1.15))
    ax.legend(fontsize=8);ax.grid(alpha=.18)
    fig.text(.5,.005,"Müller/Yan/Raoof/Weigand (2025), DaRUS-4683, CC BY 4.0. Optical response, not a volume calibration.",ha="center",fontsize=8)
    fig.savefig(OUT/"measured_optical_filling.png",dpi=160)
    fig.savefig(OUT/"measured_optical_filling.svg")
    plt.close(fig)
    out={"source":"doi:10.18419/DARUS-4683","file":chosen["filename"],
         "pixel_scale_um":2.4859,"frame_interval_s":dt,"frames":len(y),
         "nominal_pore_volume_uL":VPORE,"reported_Q1_pump_uL_s":Q,
         "result_sensitivity":results,"caveat":"Optical changes are not independently calibrated liquid volumes. No pressure synchronization."}
    (OUT/"optical_filling_results.json").write_text(json.dumps(out,indent=2))
    print("OPTICAL_FILLING_RESULTS",json.dumps(out,indent=2),flush=True)
    return out
