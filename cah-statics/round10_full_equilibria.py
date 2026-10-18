"""Round 10: physically admissible two-front saddles and secondary fold.

Planar capillary Young-Laplace model, with freely relaxing cap angles,
fixed total volume, equal pressure, and the FULL 4-DOF constrained
Hessian (2 contact positions + 2 circular-cap angles).

Do not confuse extrema on the artificially restricted Young-angle-only
fixed-volume path with stationary solutions of the full capillary problem.
All results are model predictions, not experiments.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.optimize import brentq,root
from model import PlugModel
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)
WIDTH=1.;PERIOD=.3;AMP=.05;BH=.0004;DEPTH=.38
SURFACE_TENSION=.072 # Illustrative N/m, not independently measured

class PhasePlug(PlugModel):
    def __init__(self,phase,theta_deg=36.5):
        self.phase=phase
        self.b=phase*BH
        super().__init__(theta_deg=theta_deg,width_mm=WIDTH,
                         amplitude_mm=AMP,period_mm=PERIOD,depth_mm=DEPTH)
    def w(self,x):
        return super().w(x)+self.b*(np.cos(2*self.k*np.asarray(x))-
                                     np.cos(4*self.k*np.asarray(x)))
    def dw(self,x):
        return super().dw(x)-self.b*self.k*(
            2*np.sin(2*self.k*np.asarray(x))-
            4*np.sin(4*self.k*np.asarray(x)))
    def ddw(self,x):
        return super().ddw(x)-self.b*self.k**2*(
            4*np.cos(2*self.k*np.asarray(x))-
            16*np.cos(4*self.k*np.asarray(x)))
    def B(self,x):
        x=np.asarray(x)
        return super().B(x)+self.b*(
            np.sin(2*self.k*x)/(2*self.k)-
            np.sin(4*self.k*x)/(4*self.k))
    def H(self,x):
        return self.B(x)+self.w(x)**2*self.f(self.th+self.alpha(x))
    def Hp(self,x):
        w=self.w(x);dw=self.dw(x);t=self.th+self.alpha(x)
        return (w+2*w*dw*self.f(t)+w*w*self.fp(t)*
                self.ddw(x)/(2*(1+(dw/2)**2)))
    def critical(self):
        s=brentq(lambda x:self.pprime(x,"right"),.675,.725,xtol=1e-14)
        return float(s),float(2*self.depth*self.H(s))
    def natural_state(self,c,s):
        l=c-s;r=c+s
        return np.array([l,r,self.theta(l,"left"),self.theta(r,"right")],dtype=float)
    def station(self,V,seed):
        target=V/self.depth
        def equations(q):
            c,s=q
            return [self.H(s-c)+self.H(s+c)-target,
                    self.p(s-c,"right")-self.p(s+c,"right")]
        result=root(equations,seed,tol=1e-12)
        c,s=result.x
        if (np.linalg.norm(equations([c,s]))>2e-10 or not
            (-1e-7<=c<=.18 and .5<s<.8)):
            raise RuntimeError("No physically selected pressure-balanced root")
        z=self.natural_state(c,s)
        eig,res=self.constrained_stability(z)
        assert res<2e-6 and abs(self.energy_area(z)[1]-target)<2e-8
        return dict(c_um=1000*float(c),s_mm=float(s),
            energy_mm=self.energy_area(z)[0],
            eigenvalues=[float(x) for x in eig],pressure_residual=float(
                abs(self.p(s-c,"right")-self.p(s+c,"right"))),
            kkt_residual=float(res))
    def symmetric(self,V,s0):
        s=brentq(lambda x:2*self.depth*self.H(x)-V,
                 s0-.025,s0+.025,xtol=1e-14)
        z=self.natural_state(0,s)
        eig,res=self.constrained_stability(z)
        assert res<2e-6
        return dict(c_um=0.,s_mm=float(s),
                    energy_mm=self.energy_area(z)[0],
                    eigenvalues=[float(e) for e in eig])
    def secondary_fold(self,seed=(.0236,.695)):
        def equations(q):
            c,s=q
            hl=self.Hp(s-c);hr=self.Hp(s+c)
            sdc=-(hr-hl)/(hr+hl)
            dp=self.p(s-c,"right")-self.p(s+c,"right")
            dpc=(self.pprime(s-c,"right")*(sdc-1)-
                 self.pprime(s+c,"right")*(sdc+1))
            return [dp,dpc]
        q=root(equations,seed,tol=1e-12).x
        c,s=q
        assert .005<c<.05 and np.linalg.norm(equations(q))<1e-8
        V=float(self.depth*(self.H(s-c)+self.H(s+c)))
        z=self.station(V,(c,s))
        assert abs(z["eigenvalues"][0])<2e-5
        return V,z
    def caps_clear_walls(self,state):
        c=state["c_um"]*.001;s=state["s_mm"]
        for xpos,side,sign in ((c+s,"right",1),(c-s,"left",-1)):
            width=float(self.w(xpos));ang=float(self.theta(xpos,side))
            if not 0<ang<np.pi/2:return False
            radius=width/(2*np.cos(ang))
            yy=np.linspace(0,width/2-2e-5,800)
            xx=xpos+sign*(radius*np.sin(ang)-
                           np.sqrt(np.maximum(radius*radius-yy*yy,0)))
            if float(np.min(self.w(xx)/2-yy)) < -1e-8:return False
        return True

def main():
    sens=[]
    for angle in (35.5,36.,36.5,37.,37.5,38.):
        minus=PhasePlug(-1,angle);plus=PhasePlug(1,angle)
        sA,V_A=minus.critical()
        sB,V_B=plus.critical()
        V_fold,z_fold=plus.secondary_fold()
        assert 0<V_B-V_fold<.006
        sens.append(dict(angle_deg=angle,A_onset_uL=V_A,
            B_onset_uL=V_B,B_fold_uL=V_fold,
            B_fold_c_um=z_fold["c_um"],
            B_continuous_volume_interval_nL=1000*(V_B-V_fold),
            B_fold_minus_A_onset_nL=1000*(V_fold-V_A)))
    with (OUT/"round10_wetting_sensitivity.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(sens[0]))
        writer.writeheader();writer.writerows(sens)
    mA=PhasePlug(-1);mB=PhasePlug(1)
    sA,VA=mA.critical();sB,VB=mB.critical()
    VF,fold=mB.secondary_fold()
    central_A=mA.symmetric(VA,sA);central_B=mB.symmetric(VB,sB)
    far_A=mA.station(VA,(.15,.59));far_B=mB.station(VB,(.15,.59))
    saddle_A=mA.station(VA+.001,(.011,.699))
    central_A_pre=mA.symmetric(VA+.001,sA)
    saddle_B=mB.station(VB,(.034,.692))
    emerge_B=mB.station(VB-.001,(.011,.697))
    assert abs(saddle_A["c_um"]-11.269843)<.03
    assert abs(saddle_B["c_um"]-33.629831)<.03
    assert len([x for x in saddle_A["eigenvalues"] if x< -1e-5])==1
    assert len([x for x in saddle_B["eigenvalues"] if x< -1e-5])==1
    assert all(x>0 for x in emerge_B["eigenvalues"])
    assert all(x>0 for x in far_A["eigenvalues"]+far_B["eigenvalues"])
    for m,states in ((mA,(central_A,far_A,saddle_A)),
                     (mB,(central_B,far_B,saddle_B,emerge_B,fold))):
        assert all(m.caps_clear_walls(z) for z in states)
    # Explicit counterexample: zero derivative of the artificially
    # Young-angle-restricted E(c) does NOT imply shared curvature.
    c=.03341699736899
    s=brentq(lambda x:mA.H(x-c)+mA.H(x+c)-VA/mA.depth,.67,.70)
    Hprime=float(mA.Hp(s-c))
    pdiff=float(mA.p(s-c,"right")-mA.p(s+c,"right"))
    assert abs(Hprime)<3e-6 and abs(pdiff)>.05
    pJ_per_model_mm=SURFACE_TENSION*DEPTH*1e6
    EbarA=(saddle_A["energy_mm"]-central_A_pre["energy_mm"])*pJ_per_model_mm
    EbarB=(saddle_B["energy_mm"]-central_B["energy_mm"])*pJ_per_model_mm
    assert .40<EbarA<.46 and 6.8<EbarB<7.2
    output=dict(status="MODEL ONLY; not measured or a selected dynamic landing",
       parameters=dict(width_mm=WIDTH,period_mm=PERIOD,
         main_amplitude_mm=AMP,higher_harmonic_um=.4,
         sidewall_Young_angle_deg=36.5,depth_mm=DEPTH,
         illustrative_surface_tension_N_per_m=SURFACE_TENSION),
       A=dict(onset_volume_uL=VA,at_onset=central_A,
          distant_stable=far_A,
          saddle_1_nL_above_onset=saddle_A,
          barrier_1_nL_above_onset_pJ=float(EbarA)),
       B=dict(onset_volume_uL=VB,at_onset=central_B,
          index_one_saddle_at_onset=saddle_B,
          barrier_at_onset_pJ=float(EbarB),
          stable_asymmetric_1_nL_below_onset=emerge_B,
          secondary_fold_volume_uL=VF,
          secondary_fold=fold,distant_stable=far_B),
       false_restricted_minimum=dict(centre_um=1000*c,
          pressure_mismatch_inverse_mm=pdiff,
          Hprime_one_contact=Hprime,
          explanation="E'_Young=0 due to H'=0, but p_left != p_right"),
       wetting_sensitivity=sens)
    (OUT/"round10_full_equilibria.json").write_text(json.dumps(output,indent=2))
    branches=[]
    for sign,model in ((-1,mA),(1,mB)):
        Vc=model.critical()[1]
        for c in np.linspace(.002,.041,85):
            def q(s):
                return (model.p(s-c,"right")-model.p(s+c,"right"))/c
            grid=np.linspace(.64,.725,220);yy=q(grid)
            ix=np.flatnonzero(yy[:-1]*yy[1:]<0)
            for i in ix:
                s=brentq(q,grid[i],grid[i+1])
                V=model.depth*(model.H(s-c)+model.H(s+c))
                try:z=model.station(V,(c,s))
                except RuntimeError:continue
                branches.append(dict(phase=sign,c_um=1000*c,
                    volume_offset_nL=1000*(V-Vc),
                    min_eigenvalue=z["eigenvalues"][0],
                    stable=int(z["eigenvalues"][0]>0)))
    with (OUT/"round10_bifurcation_branches.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(branches[0]))
        writer.writeheader();writer.writerows(branches)
    fig,ax=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for a,sign,name,color in zip(ax,(-1,1),("A (−0.4 µm)","B (+0.4 µm)"),
                                  ("#c57648","#0c8196")):
        for stable,style in ((1,"-"),(0,"--")):
            subset=sorted((r for r in branches if r["phase"]==sign and
                            r["stable"]==stable),key=lambda z:z["c_um"])
            if subset:a.plot([r["volume_offset_nL"] for r in subset],
                 [r["c_um"] for r in subset],style,color=color,lw=2.5,
                 label=("stable" if stable else "unstable")+" asymmetric")
        a.axvline(0,ls=":",color="0.5")
        if sign==1:a.scatter([1000*(VF-VB)],[fold["c_um"]],
                              marker="D",s=80,color=color,label="second fold")
        a.set(xlim=(-5,5),ylim=(0,48),
              xlabel="Volume relative to symmetric onset (nL)",
              ylabel="Plug-centre displacement (µm)",title="Profile "+name)
        a.grid(alpha=.18);a.legend(frameon=False,fontsize=8)
    fig.suptitle("Phase-sensitive local bifurcation and a second asymmetric fold",
                 fontweight="bold",fontsize=15)
    fig.savefig(OUT/"round10_verified_branches.png",dpi=175)
    fig.savefig(OUT/"round10_verified_branches.svg")
    print("PASS full 4-DOF constrained Hessians, equal pressure, true saddles")
    print("PASS finite-volume second asymmetric fold",VF,fold["c_um"])
    print("PASS false stationary state rejected (pressure difference)",pdiff)
    print("PASS barriers pJ, illustrative gamma",EbarA,EbarB)
    print("MODEL ONLY. No observed jumps, nonneutral 3D quartic or unique landing.")

if __name__=="__main__":
    main()
