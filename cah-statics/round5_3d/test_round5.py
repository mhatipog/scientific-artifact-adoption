"""Fast numerical controls for the finite-depth Young-Laplace study.

Checks model consistency, not experimental agreement and not continuum
stability. All test inputs are explicitly assumed ideal geometries.
"""
import numpy as np
from surface_3d import solve, K, D
from full_stability_3d import stability


def run():
    ths=np.deg2rad(36.5)
    for plate in (90,85,80):
        r=solve(.7,plate_deg=plate,a=0,b=0,deg_e=3,deg_z=2,ne=20,nz=8)
        expected=2*np.cos(ths)+2*np.cos(np.deg2rad(plate))/D
        assert r['opt'].success and r['error']<1e-4
        np.testing.assert_allclose(r['p'],expected,atol=3e-8)
    for b in (-.0005,.0005):
        r=solve(.7,plate_deg=90,b=b,deg_e=3,deg_z=2,ne=28,nz=10)
        w=1+.05*np.cos(K*.7)+b*(np.cos(2*K*.7)-np.cos(4*K*.7))
        dw=-.05*K*np.sin(K*.7)+b*(-2*K*np.sin(2*K*.7)+4*K*np.sin(4*K*.7))
        expected=2*np.cos(ths+np.arctan(dw/2))/w
        np.testing.assert_allclose(r['p'],expected,atol=5e-8)
        q=solve(.7,plate_deg=80,b=b,deg_e=3,deg_z=2,ne=28,nz=10)
        assert q['p']>r['p'] and q['error']<2e-5
        eig=stability(q,ne=28,nz=10,max_e=6,max_z=4)
        assert eig['min_total']>0 and eig['max_general_residual']<2e-4
        print(f"PASS b={b} plate80 p={q['p']:.9f} all-parity Hessian positive in tested Galerkin space")
    print("PASS independent straight-channel and planar neutral Young-Laplace limits")
    print("NOT CLAIMED: unrestricted continuum spectrum, experiment or low-angle bifurcation threshold")


if __name__=="__main__":
    run()
