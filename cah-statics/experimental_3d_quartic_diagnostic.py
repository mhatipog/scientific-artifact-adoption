"""Experimental two-front nonlinear-3D quartic diagnostic. MANUAL ONLY.

Imports the independently developed 3D variational surface energy, but
computes the quartic through a new bordered KKT implicit differentiation,
NOT a fourth-power least-squares fit. Uses reflected positive-x cap
coordinates for the left and right menisci in the even channel:
E=E_R(Y1)+E_R(Y2), V=V_R(Y1)+V_R(Y2),
C=(edge_average(Y2)-edge_average(Y1))/2.
The two faces share a total-volume constraint, not individual volumes.

A base KKT root is NOT a verified critical bifurcation. Require B2≈0,
correct 3D constrained spectrum, adaptive mesh convergence and an
independent neutral-plate continuum B4 comparison BEFORE accepting B4.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("JAX_PLATFORM_NAME","cpu")
os.environ.setdefault("XLA_FLAGS","--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
import argparse,json
import numpy as np
import jax
jax.config.update("jax_enable_x64",True)
import jax.numpy as jnp
from scipy.optimize import root
from nonlinear3d_variational import build, initial
from experimental_quartic_kkt import quartic,analytic_toy_test

S_BY_PHASE={-1:.7001915575712277,1:.6994557892661718}


def diagnostic(phase,nt,nz,plate_angle=90.):
    eps=float(np.cos(np.deg2rad(plate_angle)))
    F,dE,dV,t,z=build(phase,nt,nz)
    N=nt*nz
    s=S_BY_PHASE[phase]
    cap=initial(phase,t,z,s).ravel()
    weights=np.ones(nz)
    weights[[0,-1]]=.5
    weights/=nz-1
    weights=jnp.array(weights)
    def edge_mean(y):
        return jnp.sum((y.reshape(nt,nz)[0,:]+y.reshape(nt,nz)[-1,:])*weights)/2
    def EVC(y):
        a=y[:N];b=y[N:]
        e1,v1=F(a,eps);e2,v2=F(b,eps)
        return e1+e2,v1+v2,(edge_mean(b)-edge_mean(a))/2
    EVC=jax.jit(EVC)
    dEVC=jax.jit(jax.jacfwd(EVC))
    x0=np.r_[cap,cap]
    V0=float(EVC(jnp.array(x0))[1])
    g=np.array(dEVC(jnp.array(x0)))
    p0,mu0=np.linalg.lstsq(np.stack([g[1],g[2]],axis=1),g[0],rcond=None)[0]
    y0=np.r_[x0,p0,mu0]
    def R(y):
        x=y[:-2];p=y[-2];mu=y[-1]
        E,V,C=EVC(x)
        dE,dV,dC=dEVC(x)
        return jnp.concatenate([dE-p*dV-mu*dC,jnp.array([V-V0,C])])
    R=jax.jit(R);JR=jax.jit(jax.jacfwd(R))
    sol=root(lambda y:np.asarray(R(jnp.array(y))),y0,
             jac=lambda y:np.asarray(JR(jnp.array(y))),tol=5e-10)
    y=np.asarray(sol.x)
    resid=float(np.linalg.norm(np.asarray(R(jnp.array(y)))))
    if resid>1e-6:raise RuntimeError(f"Unable to track constrained branch: {resid}")
    formal=quartic(EVC,y,V0)
    contact_shift=float(np.mean(y[N:].reshape(nt,nz)[[0,-1],:])-s)
    result=dict(M=5,phase=phase,plate_angle_deg=plate_angle,
                transverse_nodes=nt,depth_nodes=nz,
                constrained_volume_mm3=V0,base_solver_residual=resid,
                mean_right_contact_shift_mm=contact_shift,**formal)
    result["certified_bifurcation"]=False
    result["missing_verifications"]=[
      "Find the mesh-dependent critical volume where B2=0",
      "Converge B4 with adaptive contact-line resolution and z mesh",
      "Check all orthogonal 3D constrained Hessian eigenvalues",
      "At 90 degrees recover d times independent analytical planar B4",
      "Compare branch to independent 3D FEM code"]
    print(json.dumps(result,indent=2))
    print("PROTOTYPE ONLY: even perfect stationarity does NOT certify B4 sign")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--phase",type=int,choices=[-1,1],default=1)
    parser.add_argument("--nt",type=int,default=9)
    parser.add_argument("--nz",type=int,default=3)
    parser.add_argument("--plate-angle",type=float,default=90)
    args=parser.parse_args()
    analytic_toy_test()
    diagnostic(args.phase,args.nt,args.nz,args.plate_angle)
