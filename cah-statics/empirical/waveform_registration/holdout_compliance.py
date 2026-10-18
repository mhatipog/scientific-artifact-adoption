"""Conditional 60°->120° cross-experiment compliance holdout (NOT a capillary-theory fit).

Prerequisites: bridge_camera_pressure.py generated experimental_camera_pressure_bridge.json
and manuscript/reference FIG5 original-sensor waveform, as archived in this branch.
Optical intensity is not independently calibrated to pore volume.
Publication-derived clock matching has unquantified systematic errors.
"""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).parent
OUT=HERE/"outputs"
OUT.mkdir(exist_ok=True)
BRIDGE=json.loads((OUT/"experimental_camera_pressure_bridge.json").read_text())
R={a["angle_deg"]:a for a in BRIDGE["results"]}
with (HERE.parent/"yan_video"/"observed_video_summary.csv").open() as f:
    OBS={int(x["experiment"].split("_")[0]):x for x in csv.DictReader(f)}
V_PORE=.0768 # 800x800x120 micrometre cuboid in uL
Q_PUMP=.01*1000/3600 # reported 0.01 mL/h -> uL/s

def interval(deg):
    row=OBS[deg]
    return float(row["t10_s"]),float(row["t90_s"])

def assumed_stored_volume(deg,scale=1.):
    a,b=interval(deg)
    endpoint_sat=1. if deg==60 else .92
    return .8*endpoint_sat*V_PORE*scale-Q_PUMP*(b-a)

def pressure_delta(deg,extra_shift=0.):
    csv_path=HERE/"reference"/f"fig5_{deg}deg_matched_fiso.csv"
    t_p=np.genfromtxt(csv_path,delimiter=",",skip_header=1)
    t,pressure=t_p[:,0],t_p[:,2] # original sensor, paper time axis
    a,b=interval(deg)
    shift=R[deg]["publication_time_minus_video_time_s"]+extra_shift
    return float(np.interp(b+shift,t,pressure)-np.interp(a+shift,t,pressure))

def run():
    report={"status":"CONDITIONAL experimental consistency audit; NOT direct measurement of compliance",
            "fitted_model_parameters":[],
            "model":"single shared constant compliance C (testable null hypothesis), not capillary-statics",
            "assumptions":["publication-derived camera/sensor time shifts",
                           "nominal pore cuboid volume","optical 10-to-90% corresponds to 80% nominal pore liquid volume",
                           "120 degree final saturation approximated as 92%",
                           "pressure change equals release from a shared constant apparatus compliance"],
            "experiments":{}}
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.4),layout="constrained")
    color={60:"#397896",120:"#b8754a"}
    for deg in (60,120):
        q=assumed_stored_volume(deg)
        dp=pressure_delta(deg)
        if dp>=0:raise AssertionError("Expected decreasing pressure under registered time axes")
        shifts=np.linspace(-.4,.4,161)
        drops=np.array([pressure_delta(deg,s) for s in shifts])
        valid=drops<-.01
        compliance=np.where(valid,q/abs(drops),np.nan)
        report["experiments"][str(deg)]={"nominal_storage_release_uL":q,
            "observed_pressure_change_kPa_conditional":dp,
            "inferred_C_uL_per_kPa":q/abs(dp),
            "additional_time_shift_range_s":[-.4,.4],
            "C_range_from_additional_time_shift_uL_per_kPa":[float(np.nanmin(compliance)),float(np.nanmax(compliance))]}
        axs[1].plot(shifts,compliance,color=color[deg],lw=2,label=f"{deg}°")
    v60=report["experiments"]["60"];v120=report["experiments"]["120"]
    C=v60["inferred_C_uL_per_kPa"]
    pred=-v120["nominal_storage_release_uL"]/C
    measured=v120["observed_pressure_change_kPa_conditional"]
    report["holdout"]={"60_degree_calibration_C_uL_per_kPa":C,
        "120_degree_predicted_pressure_drop_kPa":pred,
        "120_degree_observed_conditional_drop_kPa":measured,
        "error_kPa":pred-measured,
        "predicted_to_observed_drop_ratio":abs(pred/measured),
        "interpretation":"Null model disagrees under these specific assumptions; unknown volume calibration and clock errors prevent rejection of physical constant C."}
    # Explore if allowing large independent biases could rescue a shared compliance.
    scale=np.linspace(.5,1.,51)
    ratio=[]
    for a in scale:
        for b in scale:
            q60=assumed_stored_volume(60,a)
            q120=assumed_stored_volume(120,b)
            if q60>0 and q120>0:
                ratio.append((q120/abs(measured))/(q60/abs(v60["observed_pressure_change_kPa_conditional"])))
    report["volume_identifiability"]={"relative_10_90_volume_scales_tested":[.5,1.],
        "smallest_ratio_of_120_to_60_effective_C_under_independent_scalings":float(min(ratio)),
        "interpretation":"This is a sensitivity scenario, not a confidence interval."}
    axs[0].bar([0,1],[abs(v60["observed_pressure_change_kPa_conditional"]),abs(measured)],
               color=[color[60],color[120]],width=.5,label="conditional measured signal")
    axs[0].scatter([1],[abs(pred)],s=125,c="#ab3b50",marker="x",lw=3,
                   label="held-out prediction from 60° C")
    axs[0].set(xticks=[0,1],xticklabels=["60° calibration","120° holdout"],
               ylabel="Absolute pressure change (kPa)",
               title="A · A common constant compliance misses the holdout")
    axs[0].legend(frameon=False,fontsize=8)
    axs[1].set(xlabel="Additional optical time shift (seconds)",
               ylabel="Conditional C (µL/kPa)",title="B · Synchronization sensitivity",ylim=(0,1.6))
    axs[1].grid(alpha=.2);axs[1].legend(frameon=False)
    fig.suptitle("Experimental falsification gate · not a calibrated compliance measurement",weight="bold")
    fig.savefig(OUT/"conditional_compliance_holdout.png",dpi=170)
    fig.savefig(OUT/"conditional_compliance_holdout.svg")
    plt.close(fig)
    (OUT/"conditional_compliance_holdout.json").write_text(json.dumps(report,indent=2))
    assert 3.0<report["holdout"]["predicted_to_observed_drop_ratio"]<4.5
    assert v60["C_range_from_additional_time_shift_uL_per_kPa"][1]<v120["C_range_from_additional_time_shift_uL_per_kPa"][0]
    print("PASS conditional 60->120 holdout, clock sensitivity, no model pressure fitting")
    print("CAUTION: volume calibration, geometry, synchronization and capillary model predictions NOT established.")
    print(json.dumps(report["holdout"],indent=2))
if __name__=="__main__":run()
