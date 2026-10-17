"""Round 4: independent planar Young–Laplace BVP and finite-depth benchmark.

NOT a general 3D capillary solver or proof of out-of-plane stability.
All geometry in millimetres. Recreates small model-only figures/data.
"""
import csv, json
from pathlib import Path
import numpy as np, mpmath as mp
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)
mp.mp.dps=35
ANGLE=36.5
theta=np.deg2rad(ANGLE)
k=2*np.pi/.3

def coefficients(sigma, sign):
    return np.array([.05*np.exp(-.5*(k*sigma)**2),
                     sign*.0005*np.exp(-.5*(2*k*sigma)**2),
                    -sign*.0005*np.exp(-.5*(4*k*sigma)**2)])

def wall(x, coeff):
    a,b,c=coeff
    w=1+a*np.cos(k*x)+b*np.cos(2*k*x)+c*np.cos(4*k*x)
    dw=-k*(a*np.sin(k*x)+2*b*np.sin(2*k*x)+4*c*np.sin(4*k*x))
    return w,dw

def pressure(x,coeff):
    w,dw=wall(x,coeff)
    return 2*np.cos(theta+np.arctan(dw/2))/w

def bvp_test(sign):
    coeff=coefficients(0,sign)
    s=brentq(lambda x:(pressure(x+1e-6,coeff)-pressure(x-1e-6,coeff))/2e-6,
             .685,.713,xtol=1e-12)
    W,dw=wall(s,coeff)
    actual_theta=theta+np.arctan(dw/2)
    p=pressure(s,coeff)
    R=1/p
    y=np.linspace(0,W/2,120)
    exact=lambda y:s+R*np.sin(actual_theta)-np.sqrt(R*R-y*y)
    # Deliberately perturb initial shape and curvature; solve the ODE.
    initial=exact(y)+.012*np.sin(2*np.pi*y/W)**2
    ode=lambda y,u,pars:np.vstack((u[1],pars[0]*(1+u[1]**2)**1.5))
    boundary=lambda mid,end,pars:np.array([mid[1],end[0]-s,
                                           end[1]-1/np.tan(actual_theta)])
    solution=solve_bvp(ode,boundary,y,np.vstack((initial,np.gradient(initial,y))),
                       p=np.array([.94*p]),tol=4e-10,max_nodes=25000)
    assert solution.success, solution.message
    dense=np.linspace(0,W/2,2001)
    x=solution.sol(dense)[0]
    error=max(abs(x-exact(dense)))
    p_error=abs(solution.p[0]-p)
    clearance=wall(x,coeff)[0]/2-dense
    assert error<1e-8 and p_error<1e-9
    assert min(clearance[:-1])>0
    return dict(sign=sign,s=s,p=p,p_bvp=float(solution.p[0]),
                max_shape_error_mm=float(error),p_error_per_mm=float(p_error),
                minimum_sampled_interior_clearance_mm=float(min(clearance[:-1])),
                y=dense,x=x)

def quartic(sigma_mm,sign):
    sigma=mp.mpf(str(sigma_mm))
    kk=2*mp.pi/mp.mpf(".3")
    th=mp.pi*mp.mpf(str(ANGLE))/180
    a=mp.mpf(".05")*mp.exp(-mp.mpf(".5")*(kk*sigma)**2)
    b=sign*mp.mpf(".0005")*mp.exp(-mp.mpf(".5")*(2*kk*sigma)**2)
    c=-sign*mp.mpf(".0005")*mp.exp(-mp.mpf(".5")*(4*kk*sigma)**2)
    def w(x):return 1+a*mp.cos(kk*x)+b*mp.cos(2*kk*x)+c*mp.cos(4*kk*x)
    def al(x):
        return mp.atan(-kk*(a*mp.sin(kk*x)+2*b*mp.sin(2*kk*x)+
                            4*c*mp.sin(4*kk*x))/2)
    def f(t):return (2*t-mp.pi+mp.sin(2*t))/(8*mp.cos(t)**2)
    def H(x):
        return (x+a*mp.sin(kk*x)/kk+b*mp.sin(2*kk*x)/(2*kk)+
                c*mp.sin(4*kk*x)/(4*kk)+w(x)**2*f(th+al(x)))
    def p(x):return 2*mp.cos(th+al(x))/w(x)
    s=mp.findroot(lambda x:mp.diff(p,x),(mp.mpf(".69"),mp.mpf(".705")))
    return float((mp.diff(H,s,2)*mp.diff(p,s,2)-
                  mp.diff(H,s)*mp.diff(p,s,3)/3)/4)

def main():
    a,b=bvp_test(-1),bvp_test(1)
    cutoff=brentq(lambda s:quartic(s,1),.018,.021,xtol=1e-11)
    assert abs(cutoff*1000-19.9769412)<.003
    sigma=np.array([0,.5,1,2,3,5,7,10,12,15,17,19,19.5,20,22,25,30])/1000
    neg=np.array([quartic(x,-1) for x in sigma])
    pos=np.array([quartic(x,1) for x in sigma])
    assert np.all(neg<0)
    assert np.all(pos[sigma<cutoff]>0)
    assert np.all(pos[sigma>cutoff]<0)
    with (OUT/"physical_wall_smoothing.csv").open("w",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(["physical_wall_gaussian_sigma_um","B4_phase_minus","B4_phase_plus"])
        wr.writerows(zip(sigma*1000,neg,pos))
    compact=lambda result:{key:value for key,value in result.items() if key not in ("y","x")}
    (OUT/"independent_young_laplace.json").write_text(json.dumps({
        "status":"model-only, NO experimental validation",
        "exact_3d_lift":"z-independent cylinder for 90 degree top/bottom Young angle",
        "depth_note":"exact stationary solution for any constant depth; unrestricted 3D stability NOT established",
        "physical_wall_smoothing_loss_of_opposite_signs_sigma_um":cutoff*1000,
        "planar_bvp":[compact(a),compact(b)]},indent=2))
    fig,axs=plt.subplots(2,2,figsize=(13,8),layout="constrained")
    for ax,item,color in zip(axs[0],[a,b],["#ca7545","#138c8c"]):
        ax.plot(item["x"]*1000,item["y"]*1000,color=color,lw=2.4,
                label="independently solved Young–Laplace")
        ax.set(xlabel="axial meniscus position (µm)",ylabel="transverse y (µm)",
               title="Planar equilibrium, phase "+("A" if item["sign"]<0 else "B"))
        ax.grid(alpha=.2);ax.legend(frameon=False)
    axs[1,0].plot(sigma*1000,neg,"o-",color="#ca7545",label="phase A")
    axs[1,0].plot(sigma*1000,pos,"s-",color="#138c8c",label="phase B")
    axs[1,0].axhline(0,color="0.5")
    axs[1,0].axvline(cutoff*1000,ls=":",color="#7550a1",lw=2,
                     label=f"sign contrast vanishes at {cutoff*1000:.2f} µm")
    axs[1,0].set(xlabel="physical Gaussian smoothing length σ (µm)",
                  ylabel="planar quartic coefficient B₄",
                  title="How physical manufacturing smoothing changes onset")
    axs[1,0].legend(frameon=False,fontsize=8);axs[1,0].grid(alpha=.2)
    axs[1,1].axis("off")
    axs[1,1].text(.02,.92,"Finite-depth benchmark (exact stationary)",fontsize=15,
                  weight="bold",transform=axs[1,1].transAxes)
    axs[1,1].text(.02,.80,
       "For neutral (90°) top/bottom plates, extruding the\\n"
       "planar constant-curvature meniscus exactly satisfies\\n"
       "the 3D Young–Laplace equation. Principal curvatures:\\n"
       "(κ, 0); pressure = γκ. This holds at any depth.\\n\\n"
       "NOT proven: stability to z-dependent modes,\\n"
       "or behaviour at other plate contact angles.\\n\\n"
       "Image blur does not physically smooth the walls.",
       fontsize=11.5,va="top",linespacing=1.7,
       transform=axs[1,1].transAxes)
    fig.suptitle("Exact-interface benchmark and physical shape sensitivity",
                 fontsize=19,weight="bold")
    fig.savefig(OUT/"round4_young_laplace.png",dpi=180)
    fig.savefig(OUT/"round4_young_laplace.svg")
    plt.close(fig)
    print("PASS independent planar YL BVP",
          [(q["sign"],q["max_shape_error_mm"],q["p_error_per_mm"]) for q in [a,b]])
    print("PASS physical wall-smoothing threshold",cutoff*1000,"µm")
    print("3D cylinder: exact STATIONARY solution at neutral plates, no full 3D stability claim")

if __name__=="__main__":
    main()
