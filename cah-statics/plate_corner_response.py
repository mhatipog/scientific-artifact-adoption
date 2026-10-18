"""First-order *3D* plate-angle response at the two critical capillary menisci.

Solves the Jacobi PDE about the exact neutral-plate cylindrical solution:
J u = lambda, J = -u_qq-u_zz-kappa**2,
u_q(±L) = ∓kappa tan(phi) u,  u_z(0) = -eps, u_z(d) = eps.
The physical mean-curvature correction is -lambda.
eps=cos(theta_plate) is small; the numerical shape has zero-mean gauge.
This is NOT nonlinear 3D Young-Laplace or a 3D snap-type verification.
All lengths in millimetres; pressures in inverse millimetres.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import diags, kron, eye, bmat, csc_matrix
from scipy.sparse.linalg import spsolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / "outputs"
OUT.mkdir(exist_ok=True)
D = .38
K = 2*np.pi/.3
TH = np.deg2rad(36.5)


def wall(x, sign):
    w = 1+.05*np.cos(K*x)+sign*.0005*(np.cos(2*K*x)-np.cos(4*K*x))
    dw = -K*(.05*np.sin(K*x)+sign*.0005*
             (2*np.sin(2*K*x)-4*np.sin(4*K*x)))
    return w,dw


def p(x, sign):
    w,dw = wall(x,sign)
    return 2*np.cos(TH+np.arctan(dw/2))/w


def geom(sign):
    s = brentq(lambda x:(p(x+1e-6,sign)-p(x-1e-6,sign))/(2e-6),
               .685,.713,xtol=1e-13)
    w,dw = wall(s,sign)
    kap = p(s,sign)
    phi = np.pi/2-TH-np.arctan(dw/2)
    return dict(s=s,kap=kap,L=phi/kap,beta=kap*np.tan(phi))


def minus_d2(n,h,beta=0):
    diagonal = np.full(n,2/h**2)
    diagonal[[0,-1]] += 2*beta/h
    lo = np.full(n-1,-1/h**2)
    hi = lo.copy()
    lo[-1] = -2/h**2
    hi[0] = -2/h**2
    return diags([lo,diagonal,hi],[-1,0,1],format="csc")


def solve(sign,nq,nz):
    g = geom(sign)
    kap,L,beta = g["kap"],g["L"],g["beta"]
    hq,hz = 2*L/(nq-1), D/(nz-1)
    q,z = np.linspace(-L,L,nq),np.linspace(0,D,nz)
    Q = minus_d2(nq,hq,beta)-kap**2*eye(nq,format="csc")
    Z = minus_d2(nz,hz)
    J = (kron(Q,eye(nz,format="csc"))+
         kron(eye(nq,format="csc"),Z)).tocsc()
    forcing = np.zeros((nq,nz))
    forcing[:,0] += 2/hz
    forcing[:,-1] += 2/hz
    wq,wz = np.ones(nq),np.ones(nz)
    wq[[0,-1]],wz[[0,-1]] = .5,.5
    weights = np.outer(wq,wz).ravel()
    augmented = bmat([[J,-csc_matrix(np.ones((nq*nz,1)))],
                       [csc_matrix(weights.reshape(1,-1)),csc_matrix((1,1))]],
                      format="csc")
    solution = spsolve(augmented,np.r_[forcing.ravel(),0])
    U = solution[:-1].reshape(nq,nz)
    multiplier = float(solution[-1])
    assert np.max(np.abs(J@solution[:-1]-multiplier-forcing.ravel())) < 1e-8
    assert abs(weights@solution[:-1]) < 1e-8
    assert np.min(np.cos(kap*q)) > 0
    return dict(phase=sign,nq=nq,nz=nz,
                planar_kappa_mm_inv=float(kap),
                curvature_correction_per_cos_mm_inv=-multiplier,
                zero_mean_gauge_max_normal_per_cos_mm=float(np.max(np.abs(U))),
                q=q,z=z,normal_per_cos=U,L=L)


def corner_regular(side_deg,plate_deg):
    """Sharp 90° corner, acute Young angles; tests normal compatibility."""
    side,plate=np.deg2rad([side_deg,plate_deg])
    return bool(np.cos(side)**2+np.cos(plate)**2<=1+1e-12)


def run():
    rows=[]
    for phase in [-1,1]:
        earlier=None
        for nq,nz in [(33,17),(65,33),(129,65)]:
            v=solve(phase,nq,nz)
            assert abs(v["curvature_correction_per_cos_mm_inv"]-2/D)<1e-9
            if earlier is not None:
                assert abs(v["zero_mean_gauge_max_normal_per_cos_mm"]-
                           earlier["zero_mean_gauge_max_normal_per_cos_mm"])<2e-4
            earlier=v
            rows.append({k:v[k] for k in ["phase","nq","nz",
                "planar_kappa_mm_inv","curvature_correction_per_cos_mm_inv",
                "zero_mean_gauge_max_normal_per_cos_mm"]})
    assert not corner_regular(36.5,36.5)
    assert not corner_regular(36.5,50)
    assert corner_regular(36.5,60) and corner_regular(36.5,89)

    with (OUT/"round5_plate_mesh_convergence.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    finest={j:solve(j,129,65) for j in (-1,1)}
    eps=np.cos(np.deg2rad(89))
    delta_k=2*eps/D
    summary={"data_status":"linearized 3D model only; no experiment fitted",
             "plate_angle_deg":89,
             "sidewall_intrinsic_angle_deg":36.5,
             "plate_corner_film_boundary_deg":53.5,
             "jacobi_solvability_curvature_shift_mm_inv":delta_k,
             "notes":["Not a full nonlinear Young-Laplace solution",
                      "No claim that the 3D quartic bifurcation signs persist",
                      "Normal displacements are zero-mean-gauge dependent",
                      "No measured chip topography or independently calibrated Young angles"],
             "phase_results":[]}
    for j in (-1,1):
        r=finest[j]
        summary["phase_results"].append(
            {"phase":j,
             "planar_curvature_mm_inv":r["planar_kappa_mm_inv"],
             "relative_curvature_change_percent":100*delta_k/r["planar_kappa_mm_inv"],
             "max_normal_change_zero_mean_gauge_um":
                 r["zero_mean_gauge_max_normal_per_cos_mm"]*eps*1000})
    (OUT/"round5_plate_corner_results.json").write_text(json.dumps(summary,indent=2))

    fig, axs=plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for axis,phase in zip(axs[0],[-1,1]):
        v=finest[phase]
        im=axis.imshow(v["normal_per_cos"]*eps*1000,origin="lower",
            aspect="auto",extent=[0,D*1000,-v["L"]*1000,v["L"]*1000],
            vmin=-1.25,vmax=1.25,cmap="coolwarm")
        axis.set(title="Phase "+("A" if phase<0 else "B")+
                 ": 89° plate perturbation",xlabel="Depth z (µm)",
                 ylabel="Meniscus arclength (µm)")
        fig.colorbar(im,ax=axis,shrink=.82,label="Gauge-fixed normal displacement (µm)")
    axis=axs[1,0]
    angles=np.array([90,89.5,89,88,87,85,80,70,60])
    for phase in (-1,1):
        baseline=finest[phase]["planar_kappa_mm_inv"]
        axis.plot(angles,baseline+2*np.cos(np.deg2rad(angles))/D,
                  "o-",label="Phase "+("A" if phase<0 else "B"))
    axis.axvspan(60,85,color="gray",alpha=.14)
    axis.set(xlim=(90,60),xlabel="Top/bottom contact angle (°)",
             ylabel="CMC curvature (mm⁻¹)",
             title="First-order only (shaded range extrapolated)")
    axis.legend(frameon=False)
    axis=axs[1,1]
    side=np.linspace(0,90,401)
    axis.fill_between(side,0,90-side,alpha=.22,label="Corner rivulets expected")
    axis.fill_between(side,90-side,90,alpha=.16,label="No Concus–Finn corner criterion")
    axis.plot(side,90-side,"k--",linewidth=1)
    axis.scatter([36.5,36.5],[36.5,89],s=50)
    axis.set(xlim=(0,90),ylim=(0,90),xlabel="Sidewall Young angle (°)",
             ylabel="Top/bottom plate angle (°)",
             title="Sharp right-angle corner compatibility")
    axis.legend(frameon=False,fontsize=8,loc="upper right")
    fig.suptitle("First-order 3D plate-wetting response and corner regime",
                 fontweight="bold",fontsize=18)
    fig.savefig(OUT/"round5_plate_corner.png",dpi=175)
    fig.savefig(OUT/"round5_plate_corner.svg")
    plt.close(fig)
    print("PASS 3 mesh refinements, Jacobi PDE, exact solvability shift 2/d")
    print("PASS corner-angle tests and CSV/JSON/PNG/SVG outputs")
    print("89 degree plate curvature correction:",delta_k,"mm^-1")
    print("MODEL ONLY: no finite plate-angle 3D pitchfork proof")


if __name__=="__main__":
    run()
