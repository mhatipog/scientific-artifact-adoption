"""Full 3D constant-mean-curvature variational solver (symmetry restricted).

Graph X(eta,z) over a z-dependent uniform transverse coordinate:
    Y=eta * w(X_edge(z))/2,  X_edge(z)=X(eta=1,z).
Surface area and wetting free energies include all four walls.
Symmetric in y and across channel depth, top/bottom same Young angle.
Liquid volume is exact under the no-self-intersection/wall-clearance condition.
Units mm, surface tensions factored out. No actual experiment fitted.
"""
import os
os.environ.setdefault('JAX_ENABLE_X64','True')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as jnp
import numpy as np
from numpy.polynomial.legendre import leggauss, Legendre
from scipy.optimize import minimize
from functools import lru_cache
D=.38; K=2*np.pi/.3

@lru_cache(None)
def make_model(ne=32,nz=10,deg_e=4,deg_z=3,arc_n=96):
    ee,we=leggauss(ne); zz,wz=leggauss(nz)
    ta,wa=leggauss(arc_n);ta=(ta+1)/2;wa=wa/2
    modes=[(i,j) for i in range(deg_e+1) for j in range(deg_z+1)]
    def basis(e,z):
        fs=[];fe=[];fz=[]
        for i,j in modes:
            P=Legendre.basis(2*i);Q=Legendre.basis(2*j)
            fs.append(P(e)*Q(z));fe.append(P.deriv()(e)*Q(z));fz.append(P(e)*Q.deriv()(z))
        return np.array(fs).T,np.array(fe).T,np.array(fz).T
    ee2,zz2=np.meshgrid(ee,zz,indexing='ij')
    basis0,be,bz=basis(ee2.ravel(),zz2.ravel())
    edge,_,edge_z=basis(np.ones_like(zz),zz)
    plate,_,_=basis(np.repeat(ee,2),np.tile([-1,1],ne))
    cfg={k:jnp.asarray(v) for k,v in dict(be=be,bz=bz,basis=basis0,
      edge=edge,edge_z=edge_z,plate=plate,eta=ee2.ravel(),
      weights=np.outer(we,wz).ravel(),wz=wz,we=we,ta=ta,wa=wa).items()}
    def value(c,a,b,theta_side,theta_plate,refx,refeta,refedge,refplate):
        def w(x):return 1+a*jnp.cos(K*x)+b*(jnp.cos(2*K*x)-jnp.cos(4*K*x))
        def dw(x):return -a*K*jnp.sin(K*x)+b*(-2*K*jnp.sin(2*K*x)+4*K*jnp.sin(4*K*x))
        def B(x):return x+a*jnp.sin(K*x)/K+b*(jnp.sin(2*K*x)/(2*K)-jnp.sin(4*K*x)/(4*K))
        def S(x):
            t=x[:,None]*cfg['ta'][None,:]
            return x*jnp.sum(jnp.sqrt(1+(dw(t)/2)**2)*cfg['wa'][None,:],axis=1)
        xx=refx+cfg['basis']@c
        xe=refeta+cfg['be']@c
        xz=(cfg['bz']@c)*2/D
        xedge=refedge+cfg['edge']@c
        xedge_z=(cfg['edge_z']@c)*2/D
        wedge=w(xedge);wedge_prime=dw(xedge)
        weg=jnp.tile(wedge,len(cfg['we']))  # eta-major, z minor
        wz_g=jnp.tile(wedge_prime,len(cfg['we']))
        ez=jnp.tile(xedge_z,len(cfg['we']))
        # surface metric cross product in (eta,z)
        cross_z=(xe*cfg['eta']*wz_g*ez-weg*xz)/2
        surface_density=jnp.sqrt((weg/2)**2+xe**2+cross_z**2)
        A=jnp.dot(cfg['weights'],surface_density)*D/2
        # Q(z) = full channel area up to xedge minus cap's circular-segment-like exclusion
        xgrid=xx.reshape((len(cfg['we']),len(cfg['wz'])))
        Q=B(xedge)+(wedge/2)*jnp.dot(cfg['we'],xgrid-xedge[None,:])
        V=jnp.dot(cfg['wz'],Q)*D/2
        wetted_sides=D*jnp.dot(cfg['wz'],S(xedge)) # 2 sides × D/2
        xplate=(refplate+cfg['plate']@c).reshape((len(cfg['we']),2))
        # at z=±1 top/bottom, edge x by Fourier evaluation at z ±1
        edgeplate=refedge+cfg['edge']@c
        # edgeplate is at Gauss z, not endpoints. Reconstruct endpoint via xplate at eta=+1 by basis
        plate_edge_ref=refedge[0]
        ends_basis=np.asarray([Legendre.basis(2*i)(1.)*Legendre.basis(2*j)(1.) for i,j in modes])
        plate_edge=plate_edge_ref + jnp.asarray(ends_basis)@c
        plate_q=B(plate_edge)+w(plate_edge)/2*(cfg['we']@(xplate-plate_edge))
        wetted_plates=jnp.sum(plate_q) # two physical horizontal plates
        E=A-jnp.cos(theta_side)*wetted_sides-jnp.cos(theta_plate)*wetted_plates
        return jnp.stack([E,V,A,wetted_sides,wetted_plates])
    fn=jax.jit(value);jac=jax.jit(jax.jacfwd(value,argnums=0))
    return modes,fn,jac,ee,zz

def solve(s,plate_deg=90,side_deg=36.5,a=.05,b=-.0005,deg_e=4,deg_z=3,ne=32,nz=10,seed=None):
    modes,fn,jac,ee,zz=make_model(ne,nz,deg_e,deg_z)
    w=lambda x:1+a*np.cos(K*x)+b*(np.cos(2*K*x)-np.cos(4*K*x))
    dw=lambda x:-a*K*np.sin(K*x)+b*(-2*K*np.sin(2*K*x)+4*K*np.sin(4*K*x))
    phi=np.deg2rad(side_deg)+np.arctan(dw(s)/2)
    R=w(s)/(2*np.cos(phi))
    cen=s+R*np.sin(phi)
    yy=ee*w(s)/2
    x0=cen-np.sqrt(R*R-yy*yy)
    xe=ee*(w(s)/2)**2/np.sqrt(R*R-yy*yy)
    xbase=np.repeat(x0,len(zz));xeta=np.repeat(xe,len(zz))
    refedge=np.full(len(zz),s)
    refplate=np.repeat(x0,2)
    args=(a,b,np.deg2rad(side_deg),np.deg2rad(plate_deg),
          jnp.asarray(xbase),jnp.asarray(xeta),jnp.asarray(refedge),jnp.asarray(refplate))
    cc=np.zeros(len(modes)) if seed is None else seed
    target=float(fn(np.zeros(len(modes)),*args)[1])
    fun=lambda c:float(fn(c,*args)[0]);gj=lambda c:np.asarray(jac(c,*args))[0]
    cons=lambda c:float(fn(c,*args)[1]-target);gc=lambda c:np.asarray(jac(c,*args))[1]
    opt=minimize(fun,cc,method='SLSQP',jac=gj,
                 constraints=[{'type':'eq','fun':cons,'jac':gc}],
                 options={'ftol':1e-13,'maxiter':320})
    c=opt.x;data=np.asarray(fn(c,*args));gr=np.asarray(jac(c,*args));
    p=-np.dot(gr[0],gr[1])/np.dot(gr[1],gr[1]);err=np.linalg.norm(gr[0]+p*gr[1])/np.linalg.norm(gr[1]);
    return dict(p=float(p),error=float(err),c=c,modes=modes,volume=target,energy=float(data[0]),
                opt=opt,params=(s,plate_deg,side_deg,a,b),ref=(R,cen,phi),
                fn=fn,jac=jac,args=args)

if __name__=='__main__':
 for a,b in [(0,0),(.05,-.0005),(.05,.0005)]:
  for t in [90,85,80]:
   r=solve(.7,plate_deg=t,a=a,b=b)
   w=1+a*np.cos(K*.7)+b*(np.cos(2*K*.7)-np.cos(4*K*.7));dw=-a*K*np.sin(K*.7)+b*(-2*K*np.sin(2*K*.7)+4*K*np.sin(4*K*.7))
   expect=2*np.cos(np.deg2rad(36.5)+np.arctan(dw/2))/w +2*np.cos(np.deg2rad(t))/.38
   print(a,b,t,'p',r['p'],'approx',expect,'neutral diff',r['p']-expect,'error',r['error'],'nit',r['opt'].nit,'maxdc',max(abs(r['c'])))