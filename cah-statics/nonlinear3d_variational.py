"""Fully nonlinear 3D capillary surface-graph solver. Research model, not empirical fit.

Coordinates r(t,z)=[X(t,z), t*w(X(t,z))/2, z]. Liquid volume
V=integral B(X) dt dz/2; free sidewall and plate contact lines.
B'(x)=w(x), S'(x)=sqrt(1+(w'(x)/2)^2).
Find equilibrium minimizing interfacial+solid wetting energy at fixed V.
Graph representation excludes overhangs and non-graph surfaces.
"""
import os
os.environ.setdefault("JAX_PLATFORM_NAME","cpu")
os.environ.setdefault("XLA_FLAGS","--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2")
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.fft import rfft
import jax
jax.config.update("jax_enable_x64",True)
import jax.numpy as jnp

D=.38; K=2*np.pi/.3; TH=np.deg2rad(36.5); GL=(-1/np.sqrt(3),1/np.sqrt(3))
OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)

def wall(x,sgn,straight=False):
    if straight:return jnp.ones_like(x),jnp.zeros_like(x)
    a=.05;b=sgn*.0005
    return (1+a*jnp.cos(K*x)+b*(jnp.cos(2*K*x)-jnp.cos(4*K*x)),
      -K*(a*jnp.sin(K*x)+b*(2*jnp.sin(2*K*x)-4*jnp.sin(4*K*x))))

def B(x,sgn,straight=False):
    if straight:return x
    return x+.05*jnp.sin(K*x)/K+sgn*.0005*(
        jnp.sin(2*K*x)/(2*K)-jnp.sin(4*K*x)/(4*K))

def arclength_fourier(sgn,straight):
    if straight:return 1.,np.array([]),np.array([])
    x=np.arange(4096)*.3/4096
    dw=-K*(.05*np.sin(K*x)+sgn*.0005*(2*np.sin(2*K*x)-4*np.sin(4*K*x)))
    c=rfft(np.sqrt(1+(dw/2)**2))/len(x)
    n=np.arange(1,65)
    return c[0].real,2*c[1:65].real/(n*K),2*c[1:65].imag/(n*K)

def build(sgn,nt,nz,straight=False):
    c0,cr,ci=arclength_fourier(sgn,straight)
    cr=jnp.array(cr);ci=jnp.array(ci);m=jnp.arange(1,len(cr)+1,dtype=jnp.float64)
    ts=np.linspace(-1,1,nt);zs=np.linspace(0,D,nz)
    dt=2/(nt-1);dz=D/(nz-1)
    cells=[]
    for gu in GL:
      for gv in GL:
        a=(1+gu)/2;b=(1+gv)/2
        cells.append((jnp.array([(1-a)*(1-b),a*(1-b),(1-a)*b,a*b]),
            jnp.array([-(1-b),1-b,-b,b])/dt,
            jnp.array([-(1-a),-a,1-a,a])/dz,a))
    def S(x):
        if straight:return x
        phase=m*K*x[...,None]
        return c0*x+jnp.sum(cr*jnp.sin(phase)+ci*(jnp.cos(phase)-1),axis=-1)
    def energy_and_volume(X,eps):
        corners=jnp.stack([X[:-1,:-1],X[1:,:-1],X[:-1,1:],X[1:,1:]],axis=0)
        area=0.;volume=0.
        for N,Nt,Nz,a in cells:
            x=jnp.einsum("i,ijk->jk",N,corners)
            xt=jnp.einsum("i,ijk->jk",Nt,corners)
            xz=jnp.einsum("i,ijk->jk",Nz,corners)
            w,dw=wall(x,sgn,straight)
            tau=jnp.array(ts[:-1])+dt*a
            yt=w/2+tau[:,None]*dw*xt/2
            area=area+jnp.sum(jnp.sqrt(xt**2+yt**2+(w*xz/2)**2))*dt*dz/4
            volume=volume+jnp.sum(B(x,sgn,straight))*dt*dz/8
        side=0.;plate=0.
        for g in GL:
            a=(1+g)/2
            left=(1-a)*X[0,:-1]+a*X[0,1:]
            right=(1-a)*X[-1,:-1]+a*X[-1,1:]
            side=side+jnp.sum(S(left)+S(right))*dz/2
            bottom=(1-a)*X[:-1,0]+a*X[1:,0]
            top=(1-a)*X[:-1,-1]+a*X[1:,-1]
            plate=plate+jnp.sum(B(bottom,sgn,straight)+B(top,sgn,straight))*dt/4
        return area-jnp.cos(TH)*side-eps*plate,volume
    F=jax.jit(lambda x,e:energy_and_volume(x.reshape(nt,nz),e))
    dE=jax.jit(jax.grad(lambda x,e:F(x,e)[0],argnums=0))
    dV=jax.jit(jax.grad(lambda x,e:F(x,e)[1],argnums=0))
    return F,dE,dV,ts,zs

def initial(sgn,ts,zs,s,straight=False):
    if straight:
        w=1.;angle=TH;width=lambda x:1.
    else:
        w=1+.05*np.cos(K*s)+sgn*.0005*(np.cos(2*K*s)-np.cos(4*K*s))
        dw=-K*(.05*np.sin(K*s)+sgn*.0005*(2*np.sin(2*K*s)-4*np.sin(4*K*s)))
        angle=TH+np.arctan(dw/2)
        width=lambda x:1+.05*np.cos(K*x)+sgn*.0005*(np.cos(2*K*x)-np.cos(4*K*x))
    radius=w/(2*np.cos(angle))
    vals=[]
    for t in ts:
        if abs(t)==1:vals.append(s);continue
        x=s-radius*(1-np.sin(angle))
        for _ in range(100):
            y=t*width(x)/2
            nxt=s+radius*np.sin(angle)-np.sqrt(max(radius**2-y*y,0.))
            if abs(nxt-x)<1e-13:break
            x=.55*x+.45*nxt
        vals.append(nxt)
    return np.broadcast_to(np.array(vals)[:,None],(len(ts),len(zs))).copy()

def solve(sgn,nt,nz,angles,straight=False,s=.699):
    F,dE,dV,ts,zs=build(sgn,nt,nz,straight)
    guess=initial(sgn,ts,zs,s,straight).ravel()
    target=float(F(jnp.array(guess),0.)[1])
    rows=[]
    for angle in angles:
        eps=float(np.cos(np.deg2rad(angle)))
        energy=lambda x:float(F(jnp.array(x),eps)[0])
        grad=lambda x:np.asarray(dE(jnp.array(x),eps))
        constraint=lambda x:float(F(jnp.array(x),eps)[1])-target
        jac=lambda x:np.asarray(dV(jnp.array(x),eps))
        opt=minimize(energy,guess,method="SLSQP",jac=grad,
            constraints=[dict(type="eq",fun=constraint,jac=jac)],
            options=dict(ftol=1e-12,maxiter=300))
        guess=np.asarray(opt.x)
        g=grad(guess);v=jac(guess)
        p=float(np.dot(g,v)/np.dot(v,v))
        res=float(np.linalg.norm(g-p*v)/np.sqrt(len(guess)))
        # Verify graph coordinates do not fold back on themselves.
        XX=guess.reshape(nt,nz)
        Xmid=(XX[1:,:]+XX[:-1,:])/2
        Xt=(XX[1:,:]-XX[:-1,:])/(2/(nt-1))
        tm=(ts[:-1]+ts[1:])/2
        if straight:
            wmid=np.ones_like(Xmid); dwmid=np.zeros_like(Xmid)
        else:
            wmid=1+.05*np.cos(K*Xmid)+sgn*.0005*(np.cos(2*K*Xmid)-np.cos(4*K*Xmid))
            dwmid=-K*(.05*np.sin(K*Xmid)+sgn*.0005*(
                     2*np.sin(2*K*Xmid)-4*np.sin(4*K*Xmid)))
        transverse_jacobian=wmid/2+tm[:,None]*dwmid*Xt/2
        assert np.min(transverse_jacobian)>0, "Folded transverse coordinates"
        row=dict(phase=sgn,nt=nt,nz=nz,plate_angle_deg=angle,straight=straight,
           pressure_signed_per_mm=p,volume_mm3=target,energy_normalized=energy(guess),
           volume_residual_mm3=constraint(guess),stationarity_rms=res,
           max_depth_variation_um=float(1000*np.max(np.ptp(guess.reshape(nt,nz),axis=1))),
           minimum_transverse_coordinate_jacobian=float(np.min(transverse_jacobian)))
        assert opt.success and abs(row["volume_residual_mm3"])<1e-10 and res<1e-5,row
        rows.append(row)
    return rows

def run(research=False):
    rows=solve(1,25,13,[90,89,87],straight=True,s=.7)
    for row in rows:
        expected=-(2*np.cos(TH)+2*np.cos(np.deg2rad(row["plate_angle_deg"]))/D)
        assert abs(row["pressure_signed_per_mm"]-expected)<3e-4,row
    for phase in (-1,1):
        configs=[(33,9,[90,89]),(65,9,[90,89.5,89,88])]
        if research:configs.append((97,9,[90,89]))
        for nt,nz,angles in configs:
            rows+=solve(phase,nt,nz,angles)
    diagnostics=[]
    for phase in (-1,1):
        for nt in ([33,65,97] if research else [33,65]):
            subset=[r for r in rows if r["phase"]==phase and r["nt"]==nt and not r["straight"]]
            neutral=next(r for r in subset if r["plate_angle_deg"]==90)
            for r in subset:
                if r is neutral:continue
                eps=np.cos(np.deg2rad(r["plate_angle_deg"]))
                dp=-(r["pressure_signed_per_mm"]-neutral["pressure_signed_per_mm"])
                theory=2*eps/D
                remainder=r["energy_normalized"]-neutral["energy_normalized"]+2*eps*neutral["volume_mm3"]/D
                if nt>=65:assert abs(dp-theory)/theory<.03,(dp,theory)
                assert remainder<0
                diagnostics.append(dict(phase=phase,nt=nt,nz=r["nz"],angle=r["plate_angle_deg"],
                  pressure_increment_mm_inv=dp,first_order_increment_mm_inv=theory,
                  relative_error_percent=100*(dp/theory-1),
                  fixed_volume_quadratic_coefficient=remainder/eps**2))
    with (OUT/"round6_nonlinear3d_results.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(diagnostics[0]));w.writeheader();w.writerows(diagnostics)
    (OUT/"round6_nonlinear3d.json").write_text(json.dumps(dict(
       status="model-only nonlinear 3D single-front equilibria; no 3D bifurcation/empirical fit",
       exact_straight_channel_pressure="passed",
       first_order_fixed_volume_energy_invariance="supported",
       mesh_warning="corrugated absolute pressures not yet converged",
       results=diagnostics),indent=2))
    print("PASS straight-channel exact pressure, nonlinear 3D volume/stationarity")
    print("PASS quadratic fixed-volume energy remainder, near-neutral pressure increment")
    print("NOT VERIFIED: full 3D bifurcation quartic sign or experimental jump prediction")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--research",action="store_true")
    args=parser.parse_args()
    run(args.research)
