"""Independent image-curvature versus pressure diagnostic: *provisional*.

Yan et al. 2025 experimental TIFF frames and logger traces, archived under
CC BY 4.0 at doi:10.18419/DARUS-4683. No pressure fitted to camera data.
NO VERIFIED logger-to-camera synchronization or forward model validation.
"""
import json,csv
from pathlib import Path
import numpy as np,cv2
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT=Path(__file__).parent/"outputs"/"empirical"
PX_M=2.4859e-6
GAMMA=.072
GUESS={200:(130.,195.,48.),220:(173.,204.,168.)}

def fit_arc(k,low=40,high=100,band=7):
    im_path=next((OUT/"yan_frames").glob(f"120_Q1_{k:03d}_frame_*.png"))
    im=cv2.imread(str(im_path),cv2.IMREAD_GRAYSCALE)
    if im is None:raise RuntimeError(f"Unable to read {im_path}")
    edges=cv2.Canny(cv2.GaussianBlur(im,(5,5),1.5),low,high)
    yy,xx=np.nonzero(edges)
    x0,y0,r0=GUESS[k]
    selected=(abs(np.hypot(xx-x0,yy-y0)-r0)<band)&(
        xx>100)&(xx<350)&(yy>70)&(yy<345)
    points=np.c_[xx[selected],yy[selected]].astype(float)
    if len(points)<80:raise RuntimeError("Circle edge points insufficient")
    fit=least_squares(lambda z:np.hypot(points[:,0]-z[0],points[:,1]-z[1])-z[2],
                      GUESS[k],loss="soft_l1",f_scale=1.5)
    cx,cy,r=fit.x
    return dict(frame=k,cx=float(cx),cy=float(cy),r=float(r),
                radius_um=float(r*PX_M*1e6),median_residual_px=float(np.median(abs(fit.fun))),
                edge_count=len(points),canny_low=low,canny_high=high,
                seed_band_px=band,image=im)

def main():
    central=[fit_arc(i) for i in (200,220)]
    trials=[fit_arc(k,lo,hi,band) for k in (200,220)
            for lo,hi in ((25,70),(40,100),(60,140))
            for band in (5,7,10)]
    def clean(z):return {a:b for a,b in z.items() if a!="image"}
    with (OUT/"yan_radius_fit_robustness.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(clean(trials[0])));wr.writeheader()
        wr.writerows([clean(v) for v in trials])
    r1,r2=[x["r"]*PX_M for x in central]
    pred=GAMMA*(1/r1-1/r2)/1000
    dd=json.loads((OUT/"darus_lowflow_pressure_audit.json").read_text())
    sensor=next(x for x in dd["files"] if x["name"].startswith("120_"))
    pressure=sensor["decimated_trace"]
    j=min(range(len(pressure)-1),key=lambda q:pressure[q+1][1]-pressure[q][1])
    before=max(pressure[max(0,j-3):j+1],key=lambda z:z[1])
    after=min(pressure[j+1:min(len(pressure),j+5)],key=lambda z:z[1])
    measured=before[1]-after[1]
    assert 0.35<measured<0.45
    r1range=np.array([x["r"] for x in trials if x["frame"]==200])*PX_M
    r2range=np.array([x["r"] for x in trials if x["frame"]==220])*PX_M
    sensitivity=GAMMA*(1/r1range[:,None]-1/r2range[None,:])/1000
    low,hi=float(sensitivity.min()),float(sensitivity.max())
    result=dict(status="PROVISIONAL: no verified sensor-video synchronization",
                image_frame_ids=[4900,4920],video_time_camera_nominal_s=[4900/30,4920/30],
                pressure_logger_before_s=before[0],pressure_logger_after_s=after[0],
                pressure_logger_before_kPa=before[1],pressure_logger_after_kPa=after[1],
                source_doi="10.18419/DARUS-4683",source_license="CC BY 4.0",
                pixel_scale_um=2.4859,gamma_N_per_m=GAMMA,
                fitted_radius_um=[round(x["radius_um"],3) for x in central],
                curvature_derived_delta_kPa=pred,
                fit_parameter_sensitivity_range_kPa=[low,hi],
                measured_candidate_pressure_drop_kPa=measured,
                absolute_difference_kPa=abs(pred-measured),
                relative_difference_to_measured_percent=abs(pred-measured)/measured*100,
                pressure_sensor_reproducibility_kPa=0.2,
                verification="independent curvature and pressure inputs, but events NOT verified synchronized",
                limitation="normal/vertical meniscus curvature assumed constant between frames; flow/storage may intervene",
                model_type="Young-Laplace cross-modal reconstruction, NOT wall-driven forward prediction",
                novelty="classical Young-Laplace relation, not new mechanism")
    (OUT/"yan_crossmodal_provisional.json").write_text(json.dumps(result,indent=2))
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
    fig=plt.figure(figsize=(12,7.6),layout="constrained")
    axes=fig.subplots(2,2,height_ratios=[1.1,.85])
    for i,x in enumerate(central):
        ax=axes[0,i];ax.imshow(x["image"],cmap="gray",vmin=0,vmax=255)
        ax.add_patch(Circle((x["cx"],x["cy"]),x["r"],fill=False,
                            color="#00d2b9",ls="--",lw=2.4))
        ax.set(xlim=(65,440),ylim=(370,15));ax.axis("off")
        ax.text(.02,.97,f"Frame {4700+x['frame']} | R = {x['radius_um']:.0f} µm",
                transform=ax.transAxes,va="top",color="white",weight="bold",
                bbox=dict(facecolor="#173246",edgecolor="none",pad=5))
    ax=axes[1,0]
    ax.bar([0,1],[pred,measured],color=["#118899","#d27346"],width=.58)
    ax.errorbar([0],[pred],yerr=[[pred-low],[hi-pred]],fmt="none",
                capsize=5,color="#15334a")
    ax.errorbar([1],[measured],yerr=.2,fmt="none",capsize=5,color="#15334a")
    ax.set_ylim(0,.8);ax.set_ylabel("Pressure change magnitude (kPa)")
    ax.set_xticks([0,1],["From independent image radii",
                        "From pressure sensor (candidate)"],rotation=9)
    ax.grid(axis="y",alpha=.2)
    for i,v in enumerate([pred,measured]):ax.text(i,v+.03,f"{v:.3f}",ha="center",weight="bold")
    ax=axes[1,1];ax.axis("off")
    ax.text(.02,.93,
        "$|\\Delta P|\\approx\\gamma(1/R_1-1/R_2)$\\n\\n"+
        f"Curvature-derived: {pred:.3f} kPa\\n"+
        f"Observed candidate: {measured:.3f} kPa\\n"+
        f"Discrepancy: {abs(pred-measured):.3f} kPa\\n\\n"+
        "NO VERIFIED CLOCK SYNCHRONIZATION\\n"+
        "NOT a forward wall-shape prediction\\n"+
        "Vertical curvature assumed unchanged",
        transform=ax.transAxes,va="top",fontsize=10.4,linespacing=1.52,
        bbox=dict(facecolor="#f0f4f7",edgecolor="none",boxstyle="round,pad=.7"))
    fig.suptitle("First experimental image–pressure comparison (provisional)",
                 fontsize=17,weight="bold")
    fig.savefig(OUT/"yan_crossmodal_provisional.png",dpi=180)
    fig.savefig(OUT/"yan_crossmodal_provisional.svg")
    plt.close(fig)
    print("CROSSMODAL CHECK",json.dumps(result))
    print("PASS independent radii vs logger magnitude; NOT time-aligned validation")

if __name__=="__main__":main()
