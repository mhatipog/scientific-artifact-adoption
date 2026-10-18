"""Round 8: arbitrary-order isocorrelated widths with opposite capillary B4.

Mathematical 2D Young-Laplace construction; additionally an exact
neutral-plate (90-degree) 3D cylindrical stationary/linear-mode lift.
NOT fabrication-scale predictions or a nonneutral-plate 3D quartic solve.
"""
from pathlib import Path
from collections import defaultdict
import csv,json
import numpy as np, mpmath as mp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)
mp.mp.dps=38
K=2*mp.pi/mp.mpf(".3")
KF=float(K)
A=mp.mpf(".05")
BETA=mp.mpf(".046875")


def exact_mode_test(M,N):
    assert M%2==1 and M>2*N
    modes=(-1,1,-M,M,-M-2,M+2)
    states={(0,0)}
    for degree in range(1,N+1):
        states={(t+f,r+int(abs(f)>1)) for t,r in states for f in modes}
        assert not any(total==0 and r%2 for total,r in states)
    return len(states)


def quartic(M,b,angle,extra=False):
    b=mp.mpf(str(b));th=mp.mpf(str(angle))*mp.pi/180;n=M+2
    def w(x):return 1+A*mp.cos(K*x)+b*(mp.cos(M*K*x)-mp.cos(n*K*x))
    def wd(x):return -K*(A*mp.sin(K*x)+b*(M*mp.sin(M*K*x)-n*mp.sin(n*K*x)))
    def tw(x):return th+mp.atan(wd(x)/2)
    def p(x):return 2*mp.cos(tw(x))/w(x)
    def H(x):
        t=tw(x)
        f=(2*t-mp.pi+mp.sin(2*t))/(8*mp.cos(t)**2)
        return (x+A*mp.sin(K*x)/K+
                b*(mp.sin(M*K*x)/(M*K)-mp.sin(n*K*x)/(n*K))+
                w(x)**2*f)
    def approx_dp(x):
        u=KF*x; a=float(A);b0=float(b);t0=float(th)
        width=1+a*np.cos(u)+b0*(np.cos(M*u)-np.cos(n*u))
        dw=-KF*(a*np.sin(u)+b0*(M*np.sin(M*u)-n*np.sin(n*u)))
        ddw=-KF**2*(a*np.cos(u)+b0*(M*M*np.cos(M*u)-n*n*np.cos(n*u)))
        t=t0+np.arctan(dw/2)
        return 2*(-np.sin(t)*ddw/2/(1+(dw/2)**2)*width-np.cos(t)*dw)/width**2
    grid=np.linspace(.65,.74,max(501,2*M))
    val=approx_dp(grid)
    crossings=np.flatnonzero(val[:-1]*val[1:]<0)
    assert len(crossings)==1,(M,b,len(crossings))
    i=crossings[0];sf=brentq(approx_dp,grid[i],grid[i+1],xtol=2e-14)
    s=mp.findroot(lambda x:mp.diff(p,x),(mp.mpf(str(sf))-mp.mpf("1e-6"),
                                        mp.mpf(str(sf))+mp.mpf("1e-6")))
    B4=(mp.diff(H,s,2)*mp.diff(p,s,2)-mp.diff(H,s)*mp.diff(p,s,3)/3)/4
    p2=mp.diff(p,s,2)
    if not extra:return float(s),float(B4)
    C=mp.diff(H,s)*mp.sin(tw(s))*K**4/(12*w(s)*(1+(wd(s)/2)**2))
    wp2=mp.diff(w,s,2)
    normal_amplitude=1+w(s)*wp2/(4*(1+(wd(s)/2)**2))
    return float(s),float(B4),float(C),float(p2),float(normal_amplitude)


def test_clearance(M,b,s,angle):
    k=KF;n=M+2
    def w(x):return 1+.05*np.cos(k*x)+b*(np.cos(M*k*x)-np.cos(n*k*x))
    dw=-k*(.05*np.sin(k*s)+b*(M*np.sin(M*k*s)-n*np.sin(n*k*s)))
    ta=np.deg2rad(angle)+np.arctan(dw/2)
    W=w(s);R=W/(2*np.cos(ta))
    y=np.linspace(0,W/2,2001)
    x=s+R*np.sin(ta)-np.sqrt(R**2-y**2)
    clear=w(x)/2-y
    assert np.min(clear[:-1])>0 and abs(clear[-1])<1e-11
    return float(np.min(clear[:-1]))


def main():
    theta0=brentq(lambda th:quartic(5,0,th)[1],36.5,37.5,xtol=1e-10)
    s0,z,C,p2,A0=quartic(5,0,theta0,True)
    assert abs(z)<1e-7 and p2< -100 and A0>3
    u=KF*s0
    assert abs(np.sin(2*u))>.2
    rows=[]
    for M in [15,25,51,75,101,151,201]:
        N=(M-1)//2
        if M<=51:exact_mode_test(M,min(N,9))
        b=float(BETA)/M**4
        neg=quartic(M,-b,theta0,True)
        pos=quartic(M,+b,theta0,True)
        leading=float(BETA)*C*(np.cos(M*u)-np.cos((M+2)*u))
        assert neg[1]*pos[1]<0
        assert neg[4]>0 and pos[4]>0 # nodeless planar Jacobi ground states
        # Exact cylindrical stationary surfaces with neutral horizontal plates.
        if M in (15,51,201):
            test_clearance(M,-b,neg[0],theta0)
            test_clearance(M,+b,pos[0],theta0)
        rows.append(dict(M=M,correlation_order_N=N,perturbation_nm=b*1e6,
            B4_minus=neg[1],B4_plus=pos[1],
            leading_asymptotic_plus=leading,
            plus_remainder=pos[1]-leading,
            planar_ground_mode_minus=neg[4],planar_ground_mode_plus=pos[4]))
    assert abs(rows[-1]["plus_remainder"])<15
    with (OUT/"round8_arbitrary_order.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    (OUT/"round8_arbitrary_order.json").write_text(json.dumps(dict(
        scope="arbitrary finite-order correlations and neutral-plate 3D lift ONLY",
        critical_angle_deg=theta0,critical_s_mm=s0,
        p_double_prime=p2,sensitivity_C=C,jacobi_amplitude_base=A0,
        quartic_asymptotic="± beta*C*[cos(Mu)-cos((M+2)u)]+O(1/M)",
        sample=rows),indent=2))
    fig,ax=plt.subplots(figsize=(11,5),layout="constrained")
    x=np.arange(len(rows));width=.37
    ax.bar(x-width/2,[r["B4_minus"] for r in rows],width,color="#c87749",label="Phase −")
    ax.bar(x+width/2,[r["B4_plus"] for r in rows],width,color="#0e8594",label="Phase +")
    ax.axhline(0,color=".4",lw=1)
    ax.set_xticks(x,[str(r["M"]) for r in rows])
    ax.set(xlabel="High harmonic M: matching every correlation order N < M/2",
           ylabel="B₄ (normalized model units)",
           title="Arbitrary finite-order roughness statistics do not determine local capillary onset")
    ax.legend(frameon=False)
    ax.text(.02,.02,"Mathematical continuum construction. Subatomic amplitudes at high M are not physical.",
            transform=ax.transAxes,fontsize=9,color=".4")
    fig.savefig(OUT/"round8_arbitrary_order.png",dpi=170)
    fig.savefig(OUT/"round8_arbitrary_order.svg")
    print("PASS arbitrary-N Fourier correlation lemma, calibrated capillary response")
    print("PASS neutral-plate 3D cylinder ground-mode and wall-clearance sample checks")
    print("PASS asymptotic sample M=201; theorem NOT an experimental observation")


if __name__=="__main__":
    main()
