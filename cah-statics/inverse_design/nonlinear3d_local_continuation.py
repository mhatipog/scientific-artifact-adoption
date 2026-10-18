"""Geometry-fitted finite-depth free-contact-line Young–Laplace continuation.

The exact neutral cylinder is embedded analytically at all Gauss points,
with nodes equidistributed in its physical axial coordinate. The nonlinear
volume-constrained KKT equations are solved locally, NOT by unconstrained
global Ipopt minimization. This is a graph-class model and can have
unphysical alternative branches: all output is gated for monotone t->y,
small displacement, residual, and neutral-pressure convergence.

Requires NumPy, SciPy and CasADi with its sparse Newton solver.
No network or GitHub Actions. Length units mm, gamma=1.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
import numpy as np
import casadi as ca
from scipy.fft import rfft
from scipy.optimize import brentq
from nonlinear3d_casadi_audit import D,PERIOD,K,GL,geometry_numpy,pressure_planar


def cylinder_point(t,phase,s,theta_deg,straight=False):
    """Exact 90-degree-plate cylindrical graph and its dx/dt (not mesh interpolated)."""
    w,w1=geometry_numpy(s,phase,straight)
    alpha=math.radians(theta_deg)+math.atan(float(w1)/2)
    R=float(w)/(2*math.cos(alpha))
    xc=s+R*math.sin(alpha)
    if abs(t)>=1-1e-14:
        x=s
    else:
        def F(x):
            wx,_=geometry_numpy(x,phase,straight)
            return x-xc+math.sqrt(max(R*R-(t*float(wx)/2)**2,0))
        x=brentq(F,s-R*(1-math.sin(alpha))-.05,s+.03,xtol=1e-13)
    wx,wx1=geometry_numpy(x,phase,straight)
    ww=float(wx); slope=float(wx1)
    h=math.sqrt(R*R-(t*ww/2)**2)
    der=t*ww*ww/(4*h-t*t*ww*slope)
    return x,der


def physical_arclength_nodes(nt,phase,s,theta_deg,straight=False):
    """Mirror-symmetric nodes uniform in the exact circle's axial coordinate."""
    if nt%2!=1 or nt<5:
        raise ValueError("Odd nt>=5 required")
    w0,w1=geometry_numpy(s,phase,straight)
    alpha=math.radians(theta_deg)+math.atan(float(w1)/2)
    R=float(w0)/(2*math.cos(alpha))
    center=s+R*math.sin(alpha)
    x=np.linspace(center-R,s,(nt+1)//2)
    w,_=geometry_numpy(x,phase,straight)
    t=2*np.sqrt(np.maximum(0,R*R-(x-center)**2))/w
    t[0]=0.;t[-1]=1.
    if not np.all(np.diff(t)>0):
        raise ValueError("Nonmonotone solid-wall transversality")
    return np.r_[-t[:0:-1],t]

def energy_and_grad(nt,nz,phase,theta_deg,s,straight=False):
    """Return a nonlinear NLP, energy/volume and analytic gradients."""
    if nt<3 or nz<3:raise ValueError("Mesh too small")
    ts=physical_arclength_nodes(nt,phase,s,theta_deg,straight)
    delta_t=np.diff(ts)
    dz=D/(nz-1)
    X=ca.MX.sym('dX',nt*nz)
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
    gauss_base={a:np.array([cylinder_point(float(ts[i]+delta_t[i]*a),phase,s,theta_deg,straight) for i in range(nt-1)]) for a in ((1+g)/2 for g in GL)}
    dt_cells=ca.DM(np.repeat(delta_t,nz-1))
    volume=0;area=0
    for gu in GL:
        for gv in GL:
            a=(1+gu)/2;b=(1+gv)/2
            n=((1-a)*(1-b),a*(1-b),(1-a)*b,a*b)
            ntv=(-(1-b)/dt_cells,(1-b)/dt_cells,-b/dt_cells,b/dt_cells)
            nzv=(-(1-a)/dz,-a/dz,(1-a)/dz,a/dz)
            gx,gt=gauss_base[a].T
            xpos=ca.DM(np.repeat(gx,nz-1))+n[0]*a00+n[1]*a10+n[2]*a01+n[3]*a11
            xt=ca.DM(np.repeat(gt,nz-1))+ntv[0]*a00+ntv[1]*a10+ntv[2]*a01+ntv[3]*a11
            xz=nzv[0]*a00+nzv[1]*a10+nzv[2]*a01+nzv[3]*a11
            width=w(xpos)
            yt=width/2+(ca.DM(ts[ii]+delta_t[ii]*a)*wd(xpos)*xt)/2
            area+=ca.sum1(ca.sqrt(xt**2+yt**2+(width*xz/2)**2)*dt_cells)*dz/4
            volume+=ca.sum1(B(xpos)*dt_cells)*dz/8
    side=0;plate=0
    for gauss in GL:
        a=(1+gauss)/2
        vl=s+(1-a)*X[np.arange(nz-1)]+a*X[np.arange(1,nz)]
        vr=s+(1-a)*X[(nt-1)*nz+np.arange(nz-1)]+a*X[(nt-1)*nz+np.arange(1,nz)]
        side+=ca.sum1(S(vl)+S(vr))*dz/2
        base=ca.DM(gauss_base[a][:,0])
        vb=base+(1-a)*X[np.arange(nt-1)*nz]+a*X[np.arange(1,nt)*nz]
        vt=base+(1-a)*X[np.arange(nt-1)*nz+nz-1]+a*X[np.arange(1,nt)*nz+nz-1]
        plate+=ca.sum1((B(vb)+B(vt))*ca.DM(delta_t))/4
    energy=area-math.cos(math.radians(theta_deg))*side-eps*plate
    ev_fun=ca.Function("ev",[X,eps],[energy,volume],["X","eps"],["E","V"])
    grad_fun=ca.Function("grad",[X,eps],[ca.gradient(energy,X),ca.gradient(volume,X)])
    return ev_fun,grad_fun


def geometry_gate(nodes: np.ndarray, nt: int, nz: int, phase: int,
                  s: float, theta_deg: float, straight: bool) -> dict:
    """Evaluate orientation at all cell Gauss abscissae; no folded coordinates."""
    xx=nodes.reshape(nt,nz)
    tgrid=physical_arclength_nodes(nt,phase,s,theta_deg,straight)
    min_jac=float("inf")
    for i in range(nt-1):
        dt=tgrid[i+1]-tgrid[i]
        for g in GL:
            a=(1+g)/2
            t=tgrid[i]+a*dt
            x0,xt0=cylinder_point(t,phase,s,theta_deg,straight)
            for h in GL:
                b=(1+h)/2
                shift=((1-a)*(1-b)*xx[i,:-1]+a*(1-b)*xx[i+1,:-1]
                       +(1-a)*b*xx[i,1:]+a*b*xx[i+1,1:])
                dxdt=((1-b)*(xx[i+1,:-1]-xx[i,:-1])
                       +b*(xx[i+1,1:]-xx[i,1:]))/dt
                width,slope=geometry_numpy(x0+shift,phase,straight)
                J=width/2+t*slope*(xt0+dxdt)/2
                min_jac=min(min_jac,float(np.min(J)))
    return dict(min_transverse_jacobian=min_jac,
                max_cylinder_displacement_um=1000*float(np.max(abs(nodes))),
                positive_orientation=min_jac>0)


def local_continuation(nt: int,nz: int,phase: int,s: float,
                       theta_deg: float=32.,angles=(90.,89.9),
                       straight: bool=False, displacement_gate_mm: float=.01) -> list[dict]:
    if not angles or angles[0]!=90.:
        raise ValueError("Start at neutral plate angle 90 degrees")
    ev,grad=energy_and_grad(nt,nz,phase,theta_deg,s,straight)
    ndof=nt*nz
    initial=np.zeros(ndof)
    target=float(ev(initial,0)[1])
    x=ca.MX.sym("x",ndof);pressure=ca.MX.sym("pressure")
    params=ca.MX.sym("params",2)
    ge,gv=grad(x,params[0]);energy,volume=ev(x,params[0])
    func=ca.Function("KKT",[ca.vertcat(x,pressure),params],
                     [ca.vertcat(ge-pressure*gv,volume-params[1])])
    local=ca.rootfinder("local","newton",func,
                        {"abstol":1e-11,"max_iter":40,"linear_solver":"csparse"})
    exact_p=-float(pressure_planar(s,phase,theta_deg,straight))
    guess=np.r_[initial,exact_p]
    out=[]
    for angle in angles:
        eps=math.cos(math.radians(angle))
        solution=np.asarray(local(guess,np.array([eps,target]))).ravel()
        displ=solution[:-1];p=solution[-1]
        area_E,volume_actual=ev(displ,eps)
        ge_val,gv_val=grad(displ,eps)
        ge_val=np.asarray(ge_val).ravel()
        gv_val=np.asarray(gv_val).ravel()
        rms=float(np.linalg.norm(ge_val-p*gv_val)/np.sqrt(ndof))
        geom=geometry_gate(displ,nt,nz,phase,s,theta_deg,straight)
        row=dict(phase=phase,nt=nt,nz=nz,plate_angle_deg=angle,
                 epsilon=eps,pressure_signed_per_mm=float(p),
                 exact_neutral_pressure_per_mm=exact_p,
                 energy_mm2=float(area_E),V_mm3=float(volume_actual),
                 volume_residual_mm3=float(volume_actual)-target,
                 kkt_rms=rms,**geom)
        if (not geom["positive_orientation"] or
                geom["max_cylinder_displacement_um"]>1000*displacement_gate_mm or
                rms>1e-8 or abs(row["volume_residual_mm3"])>1e-9):
            raise RuntimeError(f"Local 3D continuation failed gate: {row}")
        out.append(row)
        guess=solution
    baseline=out[0]
    for row in out[1:]:
        eps=row["epsilon"]
        row["apparent_Q_mm2"]=(row["energy_mm2"]-baseline["energy_mm2"]
                              +2*eps*target/D)/(eps**2)
        row["delta_signed_pressure_per_mm"]=(row["pressure_signed_per_mm"]
                                             -baseline["pressure_signed_per_mm"])
    return out


def constrained_hessian_min_eigenvalue(nt: int,nz: int,phase: int,s: float,
                                      theta_deg: float=32., angle_deg: float=90.,
                                      straight: bool=False) -> dict:
    """Finite P1 constrained Hessian check (NOT a continuum stability proof)."""
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import LinearOperator,eigsh
    # Solve the local stationary state from the exact neutral-cylinder seed.
    ev,gr=energy_and_grad(nt,nz,phase,theta_deg,s,straight)
    n=nt*nz
    U=ca.MX.sym("U",n);pressure=ca.MX.sym("p")
    params=ca.MX.sym("par",2)
    ge,gv=gr(U,params[0]);E,V=ev(U,params[0])
    rf=ca.rootfinder("localh","newton",
        ca.Function("KKT_h",[ca.vertcat(U,pressure),params],
            [ca.vertcat(ge-pressure*gv,V-params[1])]),
        {"abstol":1e-11,"max_iter":40,"linear_solver":"csparse"})
    target=float(ev(np.zeros(n),0)[1])
    eps=math.cos(math.radians(angle_deg))
    sol=np.asarray(rf(np.r_[np.zeros(n),-float(pressure_planar(s,phase,theta_deg,straight))],
                     np.array([eps,target]))).ravel()
    X=ca.MX.sym("X",n)
    gg,vv=gr(X,eps)
    hess=ca.Function("H",[X],[ca.jacobian(gg-float(sol[-1])*vv,X)])
    H=csc_matrix(np.asarray(hess(sol[:-1])))
    v=np.asarray(gr(sol[:-1],eps)[1]).ravel()
    v=v/np.linalg.norm(v)
    def action(x):
        projected=x-v*(v@x)
        out=H@projected
        return out-v*(v@out)+v*(v@x)
    op=LinearOperator((n,n),matvec=action)
    eigen=float(eigsh(op,k=1,which="SA",tol=1e-7,maxiter=12000,
                      return_eigenvectors=False)[0])
    return dict(phase=phase,nt=nt,nz=nz,angle_deg=angle_deg,
                min_projected_hessian_eigenvalue=eigen,
                positive_discrete_hessian=eigen>0)


def audit_rows(fast: bool=False) -> list[dict]:
    rows=[]
    for phase,s in ((-1,.7039537628915958),(1,.7088690566135175)):
        meshes=((65,5),(97,7),(129,9)) if fast else (
            (65,5),(97,7),(129,9),(193,13),(257,17))
        for nt,nz in meshes:
            a,b=local_continuation(nt,nz,phase,s)
            rows.append(dict(phase=phase,nt=nt,nz=nz,
                             p90_signed_per_mm=a["pressure_signed_per_mm"],
                             neutral_pressure_error_per_mm=a["pressure_signed_per_mm"]
                                 -a["exact_neutral_pressure_per_mm"],
                             p89_9_minus_p90=b["delta_signed_pressure_per_mm"],
                             apparent_Q_89_9_mm2=b["apparent_Q_mm2"],
                             max_displacement_um=b["max_cylinder_displacement_um"],
                             min_jacobian=b["min_transverse_jacobian"],
                             max_kkt_rms=max(a["kkt_rms"],b["kkt_rms"])))
    return rows


def independent_grid_rows() -> list[dict]:
    """Separately vary planar and depth nodes; optional expensive CPU audit."""
    rows=[]
    grids=((65,9),(129,5),(129,9),(129,17),(257,9),(257,17))
    for phase,s in ((-1,.7039537628915958),(1,.7088690566135175)):
        for nt,nz in grids:
            a,b=local_continuation(nt,nz,phase,s)
            rows.append(dict(phase=phase,nt=nt,nz=nz,
                             p90_error_per_mm=a["pressure_signed_per_mm"]
                                -a["exact_neutral_pressure_per_mm"],
                             apparent_Q_89_9_mm2=b["apparent_Q_mm2"],
                             min_transverse_jacobian=b["min_transverse_jacobian"],
                             max_kkt_rms=max(a["kkt_rms"],b["kkt_rms"])))
    return rows


def as_csv(rows: list[dict]) -> str:
    import io
    buf=io.StringIO()
    writer=csv.DictWriter(buf,fieldnames=list(rows[0]),lineterminator="\n")
    writer.writeheader();writer.writerows(rows)
    return buf.getvalue()


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--fast",action="store_true")
    parser.add_argument("--grid",action="store_true",help="Independent axial/depth refinement; expensive")
    parser.add_argument("--csv",action="store_true")
    args=parser.parse_args()
    data=independent_grid_rows() if args.grid else audit_rows(fast=args.fast)
    print(as_csv(data) if args.csv else json.dumps(data,indent=2))
