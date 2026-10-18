"""Bordered KKT implicit fourth-order coefficient, no energy-polynomial fit.

For a discretized capillary energy E(x), fixed volume V(x)=V0 and an
assigned symmetry-breaking coordinate C(x)=c, define
R=(E_x-p V_x-mu C_x,V-V0,C-c)=0.
At c=0 on a verified stationary branch, B4=mu'''(0)/24.
Implements directional JAX JVPs of R and one bordered LU factorization.
This is a numerical method, NOT a claim of mesh-converged 3D B4.
"""
import numpy as np
import jax
import jax.numpy as jnp
from scipy.linalg import lu_factor,lu_solve


def quartic(EVC,y0,volume,*,tol=1e-7):
    y=jnp.asarray(y0,dtype=jnp.float64)
    dim=len(y)-2
    def residual(q):
        x=q[:-2];p=q[-2];mu=q[-1]
        E,V,C=EVC(x)
        dE,dV,dC=jax.jacfwd(EVC)(x)
        return jnp.concatenate([dE-p*dV-mu*dC,jnp.array([V-volume,C])])
    R=jax.jit(residual)
    Jfun=jax.jit(jax.jacfwd(R))
    res=np.asarray(R(y))
    error=float(np.linalg.norm(res))
    if error>=tol:raise ValueError(f"Base KKT state not certified: residual {error}")
    jac=np.asarray(Jfun(y))
    conditioning=float(np.linalg.cond(jac))
    if conditioning>1e14 or not np.isfinite(conditioning):
        raise ValueError(f"Bordered Jacobian singular: cond={conditioning}")
    fact=lu_factor(jac)
    rhs=np.zeros(dim+2);rhs[-1]=1
    yp=lu_solve(fact,rhs)
    d=jnp.asarray(yp)
    def derivative(z,u):return jax.jvp(R,(z,),(u,))[1]
    def second(z,u,v):
        return jax.jvp(lambda x:derivative(x,u),(z,),(v,))[1]
    ypp=lu_solve(fact,-np.asarray(second(y,d,d)))
    yppj=jnp.asarray(ypp)
    third=np.asarray(jax.jvp(lambda z:second(z,d,d),(y,),(d,))[1])
    mixed=np.asarray(second(y,d,yppj))
    yppp=lu_solve(fact,-third-3*mixed)
    return dict(B2=float(yp[-1]/2),B4=float(yppp[-1]/24),
                base_residual=error,bordered_condition=conditioning,
                x_response_norm=float(np.linalg.norm(yp[:-2])))


def analytic_toy_test():
    # E=A c^4 + (k/2)(y-qc^2)^2 + (s-s0)^2; V=s; C=c.
    A=123.456;k=3.2;q=5.;s0=2.5
    def EVC(x):
        c,y,s=x
        return A*c**4+k/2*(y-q*c*c)**2+(s-s0)**2,s,c
    r=quartic(EVC,np.array([0.,0.,s0,0.,0.]),s0)
    assert abs(r["B2"])<1e-9 and abs(r["B4"]-A)<1e-7,r
    print("PASS analytic KKT/Schur quartic synthetic benchmark",r)
    return r


if __name__=="__main__":
    jax.config.update("jax_enable_x64",True)
    analytic_toy_test()
