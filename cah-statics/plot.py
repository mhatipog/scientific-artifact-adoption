"""Create research figures and JSON/CSV outputs from assumptions, without fitting data."""
from pathlib import Path
import json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import PlugModel

BASE = Path(__file__).parent
OUT = BASE/"outputs"
OUT.mkdir(exist_ok=True)

m=PlugModel()
r=m.critical_and_post()
rows=[]
for theta in range(30,61,5):
    result=PlugModel(theta).critical_and_post()
    rows.append({"assumed_young_angle_deg":theta,
                 "critical_volume_uL":result["critical_volume_uL"],
                 "predicted_length_drop_um":result["predicted_length_drop_um"],
                 "predicted_center_shift_um":result["predicted_center_shift_um"]})
with (OUT/"sensitivity.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
with (OUT/"illustrative_result.json").open("w") as f:
    json.dump({"status":"theory only, not experimentally validated",
               "geometry":"assumed sinusoid, 10% peak-to-peak variation",
               "theta_Y_deg":"assumed 40, not independently measured",
               "numerical":{k:(v.tolist() if isinstance(v,np.ndarray) else v)
                            for k,v in r.items()}}, f, indent=2)

fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout="constrained")
x=np.linspace(-.90,.90,1200)
w=m.w(x)/2
ax[0].plot(x*1000,w*1000,color="#132d4b",lw=2)
ax[0].plot(x*1000,-w*1000,color="#132d4b",lw=2)
for name,z,color,style in [("before",r["before"],"#ee7d48","--"),
                            ("after",r["after"],"#17a9aa","-")]:
    L,R=z[:2]
    ax[0].axvline(L*1000,color=color,ls=style,lw=2,label=name)
    ax[0].axvline(R*1000,color=color,ls=style,lw=2)
ax[0].set(xlabel="axial coordinate (µm)",ylabel="transverse coordinate (µm)",
          title="A | Two coupled menisci (schematic)")
ax[0].legend(frameon=False)
ax[0].text(.02,.02,"Lines mark meniscus endpoints; actual interface shape not shown",
    transform=ax[0].transAxes,fontsize=7,color="0.4")
s=np.linspace(.625,.80,240)
vol=[float(m.volume(-v,v)) for v in s]
ax[1].plot(s*1000,vol,lw=2,color="#132d4b",label="symmetric branch")
ax[1].scatter([r["half_length_at_instability_mm"]*1000],
              [r["critical_volume_uL"]],color="#ee7d48",s=70,zorder=5,
              label="full stability loss")
ax[1].axhline(r["symmetric_fold_volume_uL"],color="#a0acb4",
              ls=":",label="later volume fold")
ax[1].set(xlabel="plug half-length (µm)",ylabel="volume (µL)",
          title="B | Full stability is lost before the volume fold")
ax[1].legend(frameon=False,fontsize=8)
fig.suptitle("Model prediction, not experimental jump-length agreement",weight="bold")
fig.savefig(OUT/"two_front_mechanism.png",dpi=175)
fig.savefig(OUT/"two_front_mechanism.svg")
plt.close(fig)

fig,ax=plt.subplots(figsize=(8,4.4),layout="constrained")
ax.plot([a["assumed_young_angle_deg"] for a in rows],
        [a["predicted_length_drop_um"] for a in rows],
        marker="o",color="#17a9aa",lw=2.5)
ax.set(xlabel="assumed Young angle (degrees)",
       ylabel="predicted reduction in separation (µm)",
       title="Predicted jump sensitivity — uncalibrated sinusoidal channel")
ax.grid(alpha=.2)
fig.savefig(OUT/"angle_sensitivity.png",dpi=175)
fig.savefig(OUT/"angle_sensitivity.svg")
print("Generated model-only figures and outputs in",OUT)
