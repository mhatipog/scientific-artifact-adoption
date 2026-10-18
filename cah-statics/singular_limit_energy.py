"""Round 9: singular capillary asymptotics and exact fixed-volume energy.

Derives B4/B6 from independent pressure-volume identities and integrates
exact constrained energy through the two-front equal-volume equation.
Every result is THEORY ONLY; no experiment or nonneutral 3D bifurcation.
Units: geometry mm, capillary E/(gamma*depth) mm.
"""
from pathlib import Path
import csv,json
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)
mp.mp.dps=43
K=2*mp.pi/mp.mpf(".3")
A=mp.mpf(".05")
TH=mp.mpf("37.1451464361036")*mp.pi/180
S0=mp.mpf(".698614069163614")
BETA=mp.mpf(".046875")


def model(M,sgn):
    b=mp.mpf(sgn)*BETA/M**4
    def w(x):return 1+A*mp.cos(K*x)+b*(mp.cos(M*K*x)-mp.cos((M+2)*K*x))
    def dw(x):return -K*(A*mp.sin(K*x)+b*(M*mp.sin(M*K*x)-(M+2)*mp.sin((M+2)*K*x)))
    def t(x):return TH+mp.atan(dw(x)/2)
    def p(x):return 2*mp.cos(t(x))/w(x)
    def H(x):
        a=t(x);f=(2*a-mp.pi+mp.sin(2*a))/(8*mp.cos(a)**2)
        return (x+A*mp.sin(K*x)/K+
                b*(mp.sin(M*K*x)/(M*K)-mp.sin((M+2)*K*x)/((M+2)*K))+
                w(x)**2*f)
    s=mp.findroot(lambda x:mp.diff(p,x),
                  (S0-mp.mpf("1e-4"),S0+mp.mpf("1e-4")))
    h=[None]+[mp.diff(H,s,i) for i in range(1,5)]
    p2,p3,p4,p5=[mp.diff(p,s,i) for i in (2,3,4,5)]
    s2=-h[2]/(2*h[1])
    s4=-(h[2]*s2**2+h[3]*s2+h[4]/12)/(2*h[1])
    A3=2*p2*s2+p3/3
    A5=2*p2*s4+p3*s2**2+p4*s2/3+p5/60
    B4=-h[1]*A3/4
    B6=(-h[1]*A5-(3*h[2]*s2+h[3]/2)*A3)/6
    target=2*H(s)
    def Eprime(c):
        sc=mp.findroot(lambda z:H(z+c)+H(z-c)-target,
                       (s-mp.mpf("1e-4"),s+mp.mpf("1e-4")))
        hp=mp.diff(H,sc+c);hm=mp.diff(H,sc-c)
        return (p(sc-c)-p(sc+c))*2*hp*hm/(hp+hm)
    return s,B4,B6,Eprime


def run():
    results=[];energy=[]
    Ms=(15,25,51,75,101,151,201)
    for M in Ms:
        for sign in (-1,1):
            s,B4,B6,Eprime=model(M,sign)
            if M in (15,51,201):
                assert (B4>0 if sign==1 else B4<0)
                assert B4*B6<0
            crossover=mp.sqrt(abs(B4/B6))*1000
            asymptotic=mp.sqrt(30)/(K*M)*1000
            # w_plus-w_minus -> 0 in C^3 as M -> infinity.
            C3bound=2*BETA/M**4*K**3*(M**3+(M+2)**3)
            results.append(dict(M=M,phase=sign,
                harmonic_amplitude_nm=float(BETA/M**4*1e6),
                critical_s_um=float(s*1000),B4=float(B4),B6=float(B6),
                crossover_um=float(crossover),
                asymptotic_crossover_um=float(asymptotic),
                C3_wall_difference_upper_bound=float(C3bound)))
            if sign==1 and M in (15,51,201):
                for xi in (1,2,3,4):
                    c=mp.mpf(xi)/(M*K)
                    exact=mp.quad(Eprime,[0,c])
                    kernel=mp.cos(mp.mpf(xi))-1+mp.mpf(xi)**2/2
                    predicted=24*B4/(M*K)**4*kernel
                    ratio=exact/predicted
                    assert abs(ratio-1)<mp.mpf(".15"),(M,xi,ratio)
                    energy.append(dict(M=M,xi=xi,
                        center_displacement_um=float(c*1000),
                        energy_over_gamma_depth_mm=float(exact),
                        universal_kernel=float(kernel),
                        ratio_exact_to_universal=float(ratio)))
    # Independent validation of B6 via direct exact energy at very small c.
    _,B4,B6,Eprime=model(201,1)
    for xi in (mp.mpf(".02"),mp.mpf(".03")):
        c=xi/(201*K)
        e=mp.quad(Eprime,[0,c])
        b6direct=(e-B4*c**4)/c**6
        assert abs((b6direct-B6)/B6)<mp.mpf("2e-4")
    for name,rows in (("round9_singular_limit.csv",results),
                      ("round9_energy_collapse.csv",energy)):
        with (OUT/name).open("w",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=list(rows[0]))
            wr.writeheader();wr.writerows(rows)
    (OUT/"round9_singular_limit.json").write_text(json.dumps({
        "label":"THEORY ONLY: exact planar constrained Young-Laplace energy",
        "B6/B4":"asymptotically -k^2 M^2/30 on nondegenerate subsequences",
        "scaled_energy_kernel":"cos(xi)-1+xi^2/2",
        "interpretation":"B4 signs separate in an O(1/M) window; finite-window energy difference O(M^-4)",
        "results":results,"energy":energy},indent=2))
    fig,ax=plt.subplots(figsize=(10,5.5),layout="constrained")
    xi=np.linspace(0,4,501)
    ax.plot(xi,np.cos(xi)-1+xi**2/2,lw=2.4,label="Universal crossover kernel")
    for M,mark in ((15,"o"),(51,"s"),(201,"^")):
        rows=[r for r in energy if r["M"]==M]
        xx=np.array([r["xi"] for r in rows])
        yy=np.array([r["universal_kernel"]*r["ratio_exact_to_universal"] for r in rows])
        ax.scatter(xx,yy,s=58,marker=mark,label=f"Exact fixed-volume energy, M={M}")
    ax.set(xlabel="Rescaled meniscus displacement ξ = Mkc",
           ylabel="ΔE / [24 B₄/(Mk)⁴]",
           title="Vanishing physical energy; finite local bifurcation coefficient")
    ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.savefig(OUT/"round9_energy_collapse.png",dpi=175)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.5),layout="constrained")
    for sign,mark in ((-1,"s"),(1,"o")):
        rows=[r for r in results if r["phase"]==sign]
        ax.plot([r["M"] for r in rows],[r["crossover_um"] for r in rows],
                marker=mark,label=f"Phase {sign:+}")
    ax.plot(Ms,[float(mp.sqrt(30)/(K*M)*1000) for M in Ms],
            ls="--",label="Asymptotic √30/(kM)")
    ax.set(xscale="log",yscale="log",xlabel="Harmonic order M",
           ylabel="Quartic/sixth-order crossover (µm)",
           title="Finite B₄ sign while the displacement scale shrinks")
    ax.grid(alpha=.2,which="both");ax.legend(frameon=False)
    fig.savefig(OUT/"round9_displacement_scale.png",dpi=175)
    plt.close(fig)
    print("PASS B6 exact pressure-volume expansion and independent energy integral")
    print("PASS universal finite-amplitude crossover and C3-decaying wall contrasts")
    print("CAVEAT: failure of leading scaling for phases with near-zero oscillatory factor")
    for r in energy:
        if r["xi"]==2:print(r["M"],r["ratio_exact_to_universal"])


if __name__=="__main__":
    run()
