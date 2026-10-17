"""Unrestricted (parity-complete) Galerkin Hessian of full 3D capillary energy.

Tests perturbations odd/even across y and depth around previously optimized
reflection-symmetric capillary surface. Contact lines are free on all walls.
Caution: positivity within finite polynomial graph basis is not a full
continuum stability theorem; convergence must be tested.
"""
import os
os.environ.setdefault('JAX_ENABLE_X64','True')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as jnp
import numpy as np
from numpy.polynomial.legendre import Legendre,leggauss
from scipy.linalg import null_space,eigh
from scipy.optimize import brentq
from functools import lru_cache
from surface_3d import D,K,solve

@lru_cache(None)
def get_all(ne=36,nz=12,max_e=8,max_z=6,arc_n=96):
    ee,we=leggauss(ne); zz,wz=leggauss(nz)
    t,wga=leggauss(arc_n);at=(t+1)/2;awg=wga/2
    modes=[(i,j) for i in range(max_e+1) for j in range(max_z+1)]
    def basis(e,z):
        fs=[];fe=[];fz=[]
        for i,j in modes:
            P=Legendre.basis(i);Q=Legendre.basis(j)
            fs.append(P(e)*Q(z));fe.append(P.deriv()(e)*Q(z));fz.append(P(e)*Q.deriv()(z))
        return np.array(fs).T,np.array(fe).T,np.array(fz).T
    E,Z=np.meshgrid(ee,zz,indexing='ij');basis0,be,bz=basis(E.ravel(),Z.ravel())
    er,_,erz=basis(np.ones_like(zz),zz)
    el,_,elz=basis(-np.ones_like(zz),zz)
    ep,_,_=basis(np.repeat(ee,2),np.tile([-1.,1.],ne))
    lplate,_,_=basis(np.full(2,-1.),np.array([-1.,1.]))
    rplate,_,_=basis(np.full(2,1.),np.array([-1.,1.]))
    cfg={key:jnp.asarray(v) for key,v in dict(
       be=be,bz=bz,basis=basis0,er=er,el=el,erz=erz,elz=elz,plate=ep,
       lplate=lplate,rplate=rplate,we=we,wz=wz,
       eta=E.ravel(),weight=np.outer(we,wz).ravel(),at=at,awg=awg).items()}
    def core(c, a,b,ths,thp,refx,refe,refedge,refplate):
        def w(x):return 1+a*jnp.cos(K*x)+b*(jnp.cos(2*K*x)-jnp.cos(4*K*x))
        def wp(x):return -a*K*jnp.sin(K*x)+b*(-2*K*jnp.sin(2*K*x)+4*K*jnp.sin(4*K*x))
        def B(x):return x+a*jnp.sin(K*x)/K+b*(jnp.sin(2*K*x)/(2*K)-jnp.sin(4*K*x)/(4*K))
        def S(x):
            X=x[:,None]*cfg['at'][None,:]
            return x*jnp.sum(cfg['awg'][None,:]*jnp.sqrt(1+(wp(X)/2)**2),axis=1)
        X=refx+cfg['basis']@c
        Xe=refe+cfg['be']@c
        Xz=(cfg['bz']@c)*2/D
        right=refedge+cfg['er']@c
        left=refedge+cfg['el']@c
        Rz=(cfg['erz']@c)*2/D;Lz=(cfg['elz']@c)*2/D
        wr=w(right);wl=w(left);wpr=wp(right);wpl=wp(left)
        Yeta=(wr+wl)/4
        Yz=((cfg['eta']+1)*jnp.tile(wpr*Rz,len(cfg['we']))+
             (cfg['eta']-1)*jnp.tile(wpl*Lz,len(cfg['we'])))/4
        Ye=jnp.tile(Yeta,len(cfg['we']))
        metric=jnp.sqrt(Ye**2+Xe**2+(Xe*Yz-Ye*Xz)**2)
        area=D/2*jnp.dot(cfg['weight'],metric)
        # Green theorem area of cross-section bounded by unequal sidewall endpoints.
        Q=(B(right)+B(left)-right*wr-left*wl)/2 + Yeta*(cfg['we']@X.reshape((ne,nz)))
        volume=D/2*jnp.dot(cfg['wz'],Q)
        wet_side=D/2*jnp.dot(cfg['wz'],S(left)+S(right))
        plate_x=(refplate+cfg['plate']@c).reshape((ne,2))
        left_p=refedge[0]+cfg['lplate']@c
        right_p=refedge[0]+cfg['rplate']@c
        ye=(w(left_p)+w(right_p))/4
        Qplate=((B(right_p)+B(left_p)-right_p*w(right_p)-left_p*w(left_p))/2 +
              ye*(cfg['we']@plate_x))
        wet_plate=jnp.sum(Qplate)
        energy=area-jnp.cos(ths)*wet_side-jnp.cos(thp)*wet_plate
        return jnp.stack([energy,volume])
    fn=jax.jit(core)
    jac=jax.jit(jax.jacfwd(core,argnums=0))
    hessE=jax.jit(jax.hessian(lambda c,*a:core(c,*a)[0],argnums=0))
    hessV=jax.jit(jax.hessian(lambda c,*a:core(c,*a)[1],argnums=0))
    return modes,fn,jac,hessE,hessV,ee,zz

def stability(r,ne=36,nz=12,max_e=8,max_z=6):
    modes,fn,jac,hE,hV,ee,zz=get_all(ne,nz,max_e,max_z)
    e0,z0,ths,thp,a,b=r['params'][0],r['params'][1],r['params'][2],r['params'][1],r['params'][3],r['params'][4]
    # compute exact planar reference, same one as solver
    side_deg=r['params'][2]; plate_deg=r['params'][1];s=r['params'][0]
    w=lambda x:1+a*np.cos(K*x)+b*(np.cos(2*K*x)-np.cos(4*K*x))
    dw=lambda x:-a*K*np.sin(K*x)+b*(-2*K*np.sin(2*K*x)+4*K*np.sin(4*K*x))
    phi=np.deg2rad(side_deg)+np.arctan(dw(s)/2);rad=w(s)/(2*np.cos(phi));C=s+rad*np.sin(phi)
    y=ee*w(s)/2;x=C-np.sqrt(rad*rad-y*y)
    xe=ee*(w(s)/2)**2/np.sqrt(rad*rad-y*y)
    args=(a,b,np.deg2rad(side_deg),np.deg2rad(plate_deg),
          jnp.asarray(np.repeat(x,nz)),jnp.asarray(np.repeat(xe,nz)),
          jnp.full(nz,s),jnp.asarray(np.repeat(x,2)))
    c=np.zeros(len(modes))
    for k,(i,j) in enumerate(r['modes']):
        c[modes.index((2*i,2*j))]=r['c'][k]
    gv=np.asarray(jac(c,*args))[1]
    H=np.asarray(hE(c,*args))+r['p']*np.asarray(hV(c,*args))
    H=(H+H.T)/2
    # Positive orthonormal transformation (Hessian sign). For volume constraint,
    # exclude perturbations with linear volume derivative.
    U=null_space(gv.reshape(1,-1))
    spectrum=np.linalg.eigvalsh(U.T@H@U)
    sectors={}
    for ie in [0,1]:
      for iz in [0,1]:
        ids=np.array([k for k,(i,j) in enumerate(modes) if i%2==ie and j%2==iz])
        block=H[np.ix_(ids,ids)]
        if ie==0 and iz==0:
          uv=null_space(gv[ids].reshape(1,-1));block=uv.T@block@uv
        eigen=np.linalg.eigvalsh(block)
        sectors[f'eta_{"odd" if ie else "even"}_z_{"odd" if iz else "even"}']=float(eigen[0])
    return dict(min_total=float(spectrum[0]),sectors=sectors,
                max_general_residual=float(np.linalg.norm(np.asarray(jac(c,*args))[0]+r['p']*gv)),
                matrix_dimension=len(modes))

if __name__=='__main__':
 for b in [-.0005,.0005]:
  for p in [90,85,80]:
   for deg in [3,4]:
    result=solve(.7,plate_deg=p,b=b,deg_e=deg,deg_z=deg,ne=36,nz=12)
    data=stability(result,max_e=2*deg,max_z=2*deg)
    print('PHASE',b,'PLATE',p,'DEG',deg,'PRESS',result['p'],'RESULT',data,flush=True)