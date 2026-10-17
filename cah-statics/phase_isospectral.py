"""Round 3: isospectral wall-shape perturbations reverse local snap type.

Two channels have equal mean, min/max, RMS height/slope, period and Fourier
power spectrum. Their *relative Fourier phases* differ, so local stability
need not coincide. All dimensions are mm; volume uL; energies normalized.

Theory only. No physical experiment is fitted or validated.
"""
from pathlib import Path
import csv
import json
import numpy as np
import mpmath as mp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import PlugModel

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
mp.mp.dps = 40
THETA = 36.5
A, B, LAM = 0.05, 0.0005, 0.30


class PhaseMatchedPlug(PlugModel):
    def __init__(self, b):
        self.b = b
        super().__init__(theta_deg=THETA, width_mm=1.,
                         amplitude_mm=A, period_mm=LAM, depth_mm=.38)
    def w(self,x):
        return self.w0+self.a*np.cos(self.k*np.asarray(x))+self.b*(
            np.cos(2*self.k*np.asarray(x))-np.cos(4*self.k*np.asarray(x)))
    def dw(self,x):
        return -self.a*self.k*np.sin(self.k*np.asarray(x))+self.b*(
            -2*self.k*np.sin(2*self.k*np.asarray(x))+
             4*self.k*np.sin(4*self.k*np.asarray(x)))
    def ddw(self,x):
        return -self.a*self.k**2*np.cos(self.k*np.asarray(x))+self.b*(
            -4*self.k**2*np.cos(2*self.k*np.asarray(x))+
             16*self.k**2*np.cos(4*self.k*np.asarray(x)))
    def B(self,x):
        return self.w0*np.asarray(x)+self.a/self.k*np.sin(self.k*np.asarray(x))+self.b*(
            np.sin(2*self.k*np.asarray(x))/(2*self.k)-
            np.sin(4*self.k*np.asarray(x))/(4*self.k))


def quartic(b):
    """High precision local derivative criterion of Round 2."""
    k=2*mp.pi/mp.mpf(str(LAM))
    a=mp.mpf(str(A));b=mp.mpf(str(b));th=mp.radians(THETA)
    def w(x):return 1+a*mp.cos(k*x)+b*(mp.cos(2*k*x)-mp.cos(4*k*x))
    def alpha(x):
        dw=-a*k*mp.sin(k*x)+b*(-2*k*mp.sin(2*k*x)+4*k*mp.sin(4*k*x))
        return mp.atan(dw/2)
    def f(t):return (2*t-mp.pi+mp.sin(2*t))/(8*mp.cos(t)**2)
    def H(x):
        return (x+a*mp.sin(k*x)/k+
                b*(mp.sin(2*k*x)/(2*k)-mp.sin(4*k*x)/(4*k))+
                w(x)**2*f(th+alpha(x)))
    def p(x):return 2*mp.cos(th+alpha(x))/w(x)
    sc=mp.findroot(lambda z:mp.diff(p,z),
                   (mp.mpf(".69"),mp.mpf(".71")))
    assert mp.mpf(".68")<sc<mp.mpf(".72")
    val=(mp.diff(H,sc,2)*mp.diff(p,sc,2)-
         mp.diff(H,sc)*mp.diff(p,sc,3)/3)/4
    return float(sc),float(val)


def emerging_state(model, s0, volume, sign_b4):
    """Select local asymmetric branch by equal pressure/volume, not fitting."""
    # subcritical: unstable branch appears on high-volume side;
    # supercritical: stable branch appears on low-volume side.
    dv=1e-5 if sign_b4<0 else -1e-5
    target=volume+dv
    def s_of_c(c):
        return brentq(lambda s:model.volume(c-s,c+s)-target,
                      s0-.012,s0+.012,xtol=1e-14)
    def residual(c):
        s=s_of_c(c)
        return (model.p(c-s,"left")-model.p(c+s,"right"))/c
    grid=np.linspace(.00005,.006,170)
    yy=np.array([residual(c) for c in grid])
    changes=np.where(yy[:-1]*yy[1:]<0)[0]
    if len(changes)!=1:
        raise AssertionError(f"Expected one nearby branch, found {len(changes)}")
    idx=changes[0]
    c=brentq(residual,grid[idx],grid[idx+1],xtol=1e-14)
    s=s_of_c(c)
    state=np.array([c-s,c+s,
                    float(model.theta(c-s,"left")),
                    float(model.theta(c+s,"right"))])
    ev,res=model.constrained_stability(state)
    assert res<1e-7 and abs(model.volume(state[0],state[1])-target)<1e-10
    should_stable=sign_b4>0
    assert (min(ev)>0)==should_stable
    return dict(delta_V_nL=dv*1000,centre_um=c*1000,
                stiffness_min=float(min(ev)),stable=bool(should_stable),
                energy_normalized=float(model.energy_area(state)[0]))


def wall_admissibility(m,s):
    """Check right circular cap does not cross either sinusoidal wall."""
    w=float(m.w(s));theta=float(m.theta(s,"right"))
    radius=w/(2*np.cos(theta))
    y=np.linspace(0,w/2,2001)
    x=s+radius*np.sin(theta)-np.sqrt(np.maximum(radius**2-y**2,0))
    clear=m.w(x)/2-y
    assert min(clear[:-1])>0
    return float(min(clear[:-1]))


def run():
    models=[PhaseMatchedPlug(-B),PhaseMatchedPlug(B)]
    x=np.arange(16384)*LAM/16384
    widths=[m.w(x) for m in models]
    slopes=[m.dw(x) for m in models]
    # All conventional second-order roughness metrics agree to machine precision.
    for f in [lambda y:y.mean(),np.min,np.max,np.std]:
        np.testing.assert_allclose(f(widths[0]),f(widths[1]),rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(np.std(slopes[0]),np.std(slopes[1]),rtol=1e-12)
    power=[np.abs(np.fft.rfft(w-w.mean()))**2 for w in widths]
    np.testing.assert_allclose(power[0],power[1],rtol=1e-11,atol=1e-17)
    assert abs(min(widths[0])-.95)<1e-9 and abs(max(widths[0])-1.05)<1e-9

    rows=[]
    for m,sign in zip(models,[-1,1]):
        s,b4=quartic(sign*B)
        assert b4*sign>0, (sign,b4)  # b negative subcritical; b positive supercritical
        v=float(m.volume(-s,s))
        z=m.symmetric_state(s)
        eig,r=m.constrained_stability(z)
        assert abs(eig[0])<2e-8 and r<1e-7
        branch=emerging_state(m,s,v,b4)
        clearance=wall_admissibility(m,s)
        rows.append(dict(profile="negative phase" if sign<0 else "positive phase",
                        harmonic_amplitude_um=sign*B*1000,
                        nominal_Young_angle_deg=THETA,
                        critical_half_length_um=s*1000,
                        critical_volume_uL=v,
                        quartic_B4=b4,
                        onset_type="subcritical" if b4<0 else "supercritical",
                        near_branch_center_um=branch["centre_um"],
                        near_branch_deltaV_nL=branch["delta_V_nL"],
                        near_branch_stable=branch["stable"],
                        near_branch_min_eigenvalue=branch["stiffness_min"],
                        circular_meniscus_clearance_mm=clearance))
    with (OUT/"isospectral_phase_switch.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (OUT/"isospectral_phase_switch.json").write_text(json.dumps({
        "label":"THEORY ONLY, not experimental validation",
        "same":[ "mean width", "min/max width", "rms width", "rms wall slope",
                 "Fourier power spectrum", "period", "Young angle"],
        "different":"relative Fourier phases",
        "results":rows},indent=2))

    fig=plt.figure(figsize=(12.8,8.7),layout="constrained")
    gs=fig.add_gridspec(2,2)
    ax=fig.add_subplot(gs[0,0])
    ax.plot(x*1000,widths[0]*1000,lw=2.5,color="#dd7347",label="Phase A (−0.5 µm)")
    ax.plot(x*1000,widths[1]*1000,lw=2,ls="--",color="#118a9c",label="Phase B (+0.5 µm)")
    ax.set(title="A | Nearly indistinguishable wall widths",
           xlabel="Axial distance (µm)",ylabel="Channel width (µm)")
    ax.legend(frameon=False,fontsize=9)
    ax=fig.add_subplot(gs[0,1])
    ax.plot(x*1000,(widths[1]-widths[0])*1000,color="#554a89",lw=2.1)
    ax.axhline(0,color="0.6",lw=1)
    ax.set(title="B | Phase-dependent difference, magnified",
           xlabel="Axial distance (µm)",ylabel="Width B − width A (µm)")
    ax=fig.add_subplot(gs[1,0])
    modes=np.array([1,2,3,4])
    amplitude=np.array([A,B,0,B])*1000
    ax.bar(modes-.12,amplitude,width=.24,color="#dd7347",label="Phase A")
    ax.bar(modes+.12,amplitude,width=.24,fill=False,edgecolor="#118a9c",
           linewidth=2,hatch="//",label="Phase B")
    ax.set(xlabel="Spatial harmonic number",ylabel="Fourier magnitude (µm)",
           title="C | Identical magnitude spectra",xticks=[1,2,3,4])
    ax.legend(frameon=False,fontsize=9)
    ax=fig.add_subplot(gs[1,1])
    vals=[r["quartic_B4"] for r in rows]
    ax.bar([0,1],vals,color=["#dd7347","#118a9c"],width=.55)
    ax.axhline(0,color="0.3",lw=1)
    ax.set_xticks([0,1],["Phase A","Phase B"])
    ax.set(title="D | Opposite nonlinear stability",
           ylabel="Fourth-order energy coefficient (mm⁻³)")
    ax.text(0,vals[0]/2,"SUBCRITICAL",ha="center",color="white",weight="bold",fontsize=10)
    ax.text(1,vals[1]/2,"SUPERCRITICAL",ha="center",color="white",weight="bold",fontsize=10)
    fig.suptitle("Same roughness spectrum. Different capillary instability.",
                 fontsize=21,weight="bold")
    fig.supxlabel("MODEL prediction only · 2D circular caps · θY=36.5° · no experimental fitting",
                  fontsize=10,color="#55606e")
    fig.savefig(OUT/"isospectral_snap_design.png",dpi=185)
    fig.savefig(OUT/"isospectral_snap_design.svg")
    plt.close(fig)
    print("PASS identical widths and roughness power statistics:",[(r["profile"],r["quartic_B4"],r["near_branch_stable"]) for r in rows])
    print("PASS full constrained Hessian, adjacent pressure/volume branches, cap clearance")
    print("NOT TESTED: measured wall profiles, real pinning, jump lengths, 3D meniscus stability")


if __name__=="__main__":
    run()
