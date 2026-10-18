"""Independent CasADi finite-element free-contact-line 3D model.

Parametric surface r(t,z)=[X(t,z), t*w(X(t,z))/2, z].
Uses tensor bilinear interpolation and 2x2 Gauss area/volume quadrature.
This is a local graph-class numerical model; not a global stability proof.
"""
from __future__ import annotations
import math
import numpy as np
import casadi as ca
from scipy.fft import rfft
from scipy.optimize import brentq

D=.38
PERIOD=.3
K=2*math.pi/PERIOD
GL=(-1/math.sqrt(3),1/math.sqrt(3))

def geometry_numpy(x, phase, straight=False):
    x=np.asarray(x)
    if straight:
        return np.ones_like(x),np.zeros_like(x)
    u=K*x
    return (1+.05*np.cos(u)+phase*.00017*(np.cos(5*u)-np.cos(7*u)),
            -K*(.05*np.sin(u)+phase*.00017*(5*np.sin(5*u)-7*np.sin(7*u))))

def pressure_planar(x,phase,theta_deg,straight=False):
    w,w1=geometry_numpy(x,phase,straight)
    return 2*np.cos(np.deg2rad(theta_deg)+np.arctan(w1/2))/w

def initial_cylinder(nt,nz,phase,s,theta_deg,straight=False):
    w,w1=geometry_numpy(s,phase,straight)
    angle=math.radians(theta_deg)+math.atan(float(w1/2))
    radius=float(w)/(2*math.cos(angle))
    tvec=np.linspace(-1,1,nt)
    values=[]
    for t in tvec:
        if abs(t)==1:
            values.append(s)
            continue
        def f(x):
            wx,_=geometry_numpy(x,phase,straight)
            y=t*float(wx)/2
            return x-s-radius*math.sin(angle)+math.sqrt(max(radius*radius-y*y,0))
        low=s-radius*(1-math.sin(angle))-.05
        xx=brentq(f,low,s+.03,xtol=1e-13)
        values.append(xx)
    return np.broadcast_to(np.array(values)[:,None],(nt,nz)).copy().ravel()

def solver_model(nt,nz,phase,theta_deg,straight=False):
    """Return a nonlinear NLP, energy/volume and analytic gradients."""
    if nt<3 or nz<3:raise ValueError("Mesh too small")
    dt=2/(nt-1);dz=D/(nz-1)
    X=ca.MX.sym('X',nt*nz)
    eps=ca.MX.sym('eps')
    def w(x):
        if straight:return x*0+1
        return 1+.05*ca.cos(K*x)+phase*.00017*(ca.cos(5*K*x)-ca.cos(7*K*x))
    def wd(x):
        if straight:return x*0
        return -K*(.05*ca.sin(K*x)+phase*.00017*(5*ca.sin(5*K*x)-7*ca.sin(7*K*x)))
    def B(x):
        if straight:return x
        return x+.05*ca.sin(K*x)/K+phase*.00017*(ca.sin(5*K*x)/(5*K)-ca.sin(7*K*x)/(7*K))
    if straight:
        coeff=(1.,np.array([]),np.array([]))
    else:
        xs=np.arange(4096)*PERIOD/4096
        _,dx=geometry_numpy(xs,phase,False)
        cf=rfft(np.sqrt(1+(dx/2)**2))/len(xs)
        nums=np.arange(1,65)
        coeff=(float(cf[0].real),2*cf[1:65].real/(nums*K),
               2*cf[1:65].imag/(nums*K))
    c0,cr,ci=coeff
    def S(x):
        if straight:return x
        v=c0*x
        for ni in range(1,len(cr)+1):
            v=v+float(cr[ni-1])*ca.sin(ni*K*x)+float(ci[ni-1])*(ca.cos(ni*K*x)-1)
        return v
    ii=np.repeat(np.arange(nt-1),nz-1)
    jj=np.tile(np.arange(nz-1),nt-1)
    k00=ii*nz+jj;k10=(ii+1)*nz+jj;k01=ii*nz+jj+1;k11=(ii+1)*nz+jj+1
    a00,a10,a01,a11=[X[idx] for idx in (k00,k10,k01,k11)]
    volume=0;area=0
    for gu in GL:
        for gv in GL:
            a=(1+gu)/2;b=(1+gv)/2
            n=((1-a)*(1-b),a*(1-b),(1-a)*b,a*b)
            ntv=(-(1-b)/dt,(1-b)/dt,-b/dt,b/dt)
            nzv=(-(1-a)/dz,-a/dz,(1-a)/dz,a/dz)
            xpos=n[0]*a00+n[1]*a10+n[2]*a01+n[3]*a11
            xt=ntv[0]*a00+ntv[1]*a10+ntv[2]*a01+ntv[3]*a11
            xz=nzv[0]*a00+nzv[1]*a10+nzv[2]*a01+nzv[3]*a11
            width=w(xpos)
            yt=width/2+(ca.DM(-1+dt*(ii+a))*wd(xpos)*xt)/2
            area+=ca.sum1(ca.sqrt(xt**2+yt**2+(width*xz/2)**2))*dt*dz/4
            volume+=ca.sum1(B(xpos))*dt*dz/8
    side=0;plate=0
    for gauss in GL:
        a=(1+gauss)/2
        vl=(1-a)*X[np.arange(nz-1)]+a*X[np.arange(1,nz)]
        vr=(1-a)*X[(nt-1)*nz+np.arange(nz-1)]+a*X[(nt-1)*nz+np.arange(1,nz)]
        side+=ca.sum1(S(vl)+S(vr))*dz/2
        vb=(1-a)*X[np.arange(nt-1)*nz]+a*X[np.arange(1,nt)*nz]
        vt=(1-a)*X[np.arange(nt-1)*nz+nz-1]+a*X[np.arange(1,nt)*nz+nz-1]
        plate+=ca.sum1(B(vb)+B(vt))*dt/4
    energy=area-math.cos(math.radians(theta_deg))*side-eps*plate
    ev_fun=ca.Function("ev",[X,eps],[energy,volume],["X","eps"],["E","V"])
    grad_fun=ca.Function("grad",[X,eps],[ca.gradient(energy,X),ca.gradient(volume,X)])
    # Fixed volume enforced through lbg=ubg at solve time.
    solver=ca.nlpsol('capillary','ipopt',{'x':X,'p':eps,'f':energy,'g':volume},
                {'ipopt.print_level':0,'print_time':0,'ipopt.max_iter':450,
                 'ipopt.tol':1e-10,'ipopt.acceptable_tol':1e-9})
    return solver,ev_fun,grad_fun

def run(nt=13,nz=5,phase=-1,theta_deg=32,straight=True,angle=90.,s=.7):
    solver,ev,grad=solver_model(nt,nz,phase,theta_deg,straight)
    x0=initial_cylinder(nt,nz,phase,s,theta_deg,straight)
    e0,v0=ev(x0,0)
    eps=math.cos(math.radians(angle))
    result=solver(x0=x0,p=eps,lbg=float(v0),ubg=float(v0))
    xx=np.asarray(result['x']).ravel()
    EE,VV=ev(xx,eps)
    ge,gv=grad(xx,eps)
    gg=np.asarray(ge).ravel();vv=np.asarray(gv).ravel()
    press=float(gg@vv/(vv@vv))
    res=np.linalg.norm(gg-press*vv)/math.sqrt(len(xx))
    profile=xx.reshape(nt,nz)
    return dict(nt=nt,nz=nz,phase=phase,angle=angle,success=solver.stats()["success"],
        n_iter=solver.stats()["iter_count"],norm_E=float(EE),V=float(VV),V0=float(v0),
        pressure_signed=press,pressure_planar=pressure_planar(s,phase,theta_deg,straight),
        stationarity_rms=res,max_depth_um=1000*max(np.ptp(profile,axis=1)),
        max_init_dev_um=1000*max(abs(xx-x0)),min_X=float(min(xx)),max_X=float(max(xx)))

def solve_sequence(nt,nz,phase,theta_deg,angles,s,straight=False):
    solver,ev,grad=solver_model(nt,nz,phase,theta_deg,straight)
    initial=initial_cylinder(nt,nz,phase,s,theta_deg,straight)
    _,v0=ev(initial,0)
    target=float(v0)
    guess=initial
    results=[]
    ts=np.linspace(-1,1,nt)
    for angle in angles:
        eps=math.cos(math.radians(angle))
        ans=solver(x0=guess,p=eps,lbg=target,ubg=target)
        if not solver.stats()["success"]:
            raise RuntimeError(f"IPOPT failed {nt},{nz},{phase},{angle}")
        xx=np.asarray(ans["x"]).ravel()
        energy,vol=ev(xx,eps)
        gr,gv=grad(xx,eps)
        gr=np.asarray(gr).ravel();gv=np.asarray(gv).ravel()
        p=float(gr@gv/(gv@gv))
        rms=float(np.linalg.norm(gr-p*gv)/math.sqrt(len(xx)))
        XX=xx.reshape(nt,nz)
        Xt=(XX[1:,:]+XX[:-1,:])/2
        dXt=(XX[1:,:]-XX[:-1,:])/(2/(nt-1))
        width,wd=geometry_numpy(Xt,phase,straight)
        jac=width/2+((ts[:-1]+ts[1:])[:,None]/2)*wd*dXt/2
        row=dict(nt=nt,nz=nz,phase=phase,angle_deg=angle,
            pressure_signed_per_mm=p,pressure_analytic_neutral_per_mm=-pressure_planar(s,phase,theta_deg,straight),
            V_target_mm3=target,V_res_mm3=float(vol)-target,
            energy_normalized=float(energy),stationarity_rms=rms,
            ipopt_iterations=solver.stats()["iter_count"],
            max_depth_span_um=float(1000*np.max(np.ptp(XX,axis=1))),
            min_transverse_jac=float(jac.min()),
            max_shift_from_cylinder_um=float(1000*np.max(np.abs(xx-initial))))
        if row["min_transverse_jac"]<=0:raise RuntimeError("Folded coordinates")
        results.append(row)
        guess=xx
    return results


def convergence_rows():
    """Local P1 3D pilot. The absolute pressure has FAILED a mesh gate."""
    rows=[]
    for phase,s in ((-1,.7039537628915958),(1,.7088690566135175)):
        for nt,nz in ((17,5),(33,7),(49,9),(65,9),(97,13)):
            a,b=solve_sequence(nt,nz,phase,32,[90.,89.],s)
            eps=math.cos(math.radians(89))
            q=(b["energy_normalized"]-a["energy_normalized"]
               +2*eps*a["V_target_mm3"]/D)/(eps**2)
            rows.append(dict(phase=phase,n_t=nt,n_z=nz,
                             signed_p90_per_mm=a["pressure_signed_per_mm"],
                             exact_p90_per_mm=a["pressure_analytic_neutral_per_mm"],
                             delta_p_89minus90_per_mm=(b["pressure_signed_per_mm"]
                                 -a["pressure_signed_per_mm"]),
                             apparent_Q89_mm2=q,
                             max_kkt_rms=max(a["stationarity_rms"],
                                             b["stationarity_rms"]),
                             status="REJECT_UNCONVERGED_ABSOLUTE_PRESSURE_AND_Q"))
    return rows


if __name__=="__main__":
    import argparse,csv,sys,json
    parser=argparse.ArgumentParser()
    parser.add_argument("--audit",action="store_true",
                        help="10 corrugated mesh cases; no writes")
    args=parser.parse_args()
    if args.audit:
        rows=convergence_rows()
        writer=csv.DictWriter(sys.stdout,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    else:
        print(json.dumps(run(),indent=2))
