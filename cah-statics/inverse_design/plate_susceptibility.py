"""Neutral-plate capillary plate-angle susceptibility at a stationary cylindrical interface.

At fixed liquid volume and with eps=cos(theta_plate),
 E(eps)=E(0)-2 eps V/depth + Q eps**2 + O(eps**3).
Here Q is a *second-order energy coefficient*, not a two-front quartic.
Physical units: mm, E/gamma in mm**2, Q in mm**2.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import bmat, csc_matrix, diags, kron
from scipy.sparse.linalg import spsolve

from spectral_isolation import Channel, critical_state


def even_mode(t_index: int, L: float, beta: float) -> tuple[float, float]:
    """Even Robin eigenfunction wave number and squared spatial integral."""
    if not (t_index >= 0 and beta > 0 and L > 0):
        raise ValueError("Need nonnegative even mode, positive Robin coefficient and half-span")
    a = t_index * math.pi + 1e-12
    b = (t_index + 0.5) * math.pi - 1e-12
    t = brentq(lambda v: v * math.tan(v) - beta * L, a, b,
               xtol=1e-13)
    k = t / L
    normalization = L + math.sin(2*t) / (2*k)
    integral_sq = (2*math.sin(t)/k)**2 / normalization
    return k, integral_sq


def susceptibility_exact(state: dict, depth_mm: float, even_modes: int = 80) -> dict:
    """Exact neutral-cylinder continuum Jacobi spectral expansion.

    Includes pressure/volume projection of the zero mode. Strictly Q < 0 for
    positive depth, the specified local neutral-plate Jacobi geometry and
    contact-angle sign convention.
    """
    if depth_mm <= 0 or even_modes < 2:
        raise ValueError("Need depth>0 and at least two even modes")
    L = state["L"]
    kappa = state["kappa"]
    beta = state["beta_wall"]
    if not (0 < kappa*L < math.pi/2):
        raise ValueError("Nodeless neutral critical state required")
    if abs(beta - kappa*math.tan(kappa*L)) > 2e-9:
        raise ValueError("Wall Robin coefficient does not possess the critical zero mode")
    k0, g0_sq = even_mode(0,L,beta)
    if abs(k0 - kappa) > 5e-9:
        raise ArithmeticError("Critical root mismatch")
    zero_part = -depth_mm * g0_sq / 6
    higher = 0.0
    first_excited = None
    for j in range(1,even_modes):
        kj, gj_sq = even_mode(j,L,beta)
        lam = kj**2-kappa**2
        if lam <= 0:
            raise ArithmeticError("Nonpositive excited eigenvalue")
        root = math.sqrt(lam)
        x = root*depth_mm/2
        bracket = (1/math.tanh(x)-1/x)/root
        part = -gj_sq*bracket
        if first_excited is None:
            first_excited = part
        higher += part
    # Rigorous omitted-mode bound: for j>=N>=2, k_j >= j*pi/L,
    # N_j >= L, |sin(k_j L)| <= beta/k_j and
    # sqrt(k_j*k_j-kappa*kappa) >= sqrt(3)*k_j/2.
    # Thus |term_j| <= 8*beta**2*L**4/(sqrt(3)*pi**5*j**5).
    # The integral comparison sum_{j>=N} j**-5 <= 1/(4*(N-1)**4).
    tail_bound=2*beta**2*L**4/(math.sqrt(3)*math.pi**5*(even_modes-1)**4)
    return dict(Q_mm2=zero_part+higher,
                Q_lower_bound_mm2=zero_part+higher-tail_bound,
                Q_upper_bound_mm2=zero_part+higher,
                spectral_tail_bound_mm2=tail_bound,
                zero_mode_mm2=zero_part,
                first_excited_mm2=first_excited,
                higher_modes_mm2=higher,even_modes=even_modes,
                pressure_response_per_mm=2/depth_mm,
                initial_slope_energy_per_mm3=-2/depth_mm)


def p1_1d(nodes: int, half_length: float, endpoint_robin: float = 0):
    if nodes < 3 or half_length <= 0:
        raise ValueError("Invalid P1 domain")
    h = 2*half_length/(nodes-1)
    mdiag = np.full(nodes,2*h/3)
    mdiag[[0,-1]]=h/3
    M = diags([np.full(nodes-1,h/6), mdiag, np.full(nodes-1,h/6)],
              [-1,0,1],format="csc")
    kdiag = np.full(nodes,2/h)
    kdiag[[0,-1]]=1/h
    K = diags([np.full(nodes-1,-1/h),kdiag,np.full(nodes-1,-1/h)],
              [-1,0,1],format="csc")
    if endpoint_robin:
        K=K+diags([np.r_[endpoint_robin,np.zeros(nodes-2),endpoint_robin]],
                  [0],format="csc")
    return K,M


def susceptibility_fem(state: dict, depth_mm: float, q_nodes: int, z_nodes: int) -> dict:
    """Independent conforming tensor P1 solve with global volume constraint.

    Coordinates: q in [-L,L], z in [0,depth], with physical depth, not 2*depth.
    (The half_length helper therefore receives depth/2 for the z mesh.)
    Weak form: (K u,v)+lambda*int(v)=plate(v), int(u)=0.
    """
    if depth_mm <= 0:
        raise ValueError("Depth must be positive")
    L,kappa,beta=(state[key] for key in ("L","kappa","beta_wall"))
    Kq,Mq=p1_1d(q_nodes,L,beta)
    Kz,Mz=p1_1d(z_nodes,depth_mm/2)
    A=(kron(Kq-kappa**2*Mq,Mz)+kron(Mq,Kz)).tocsc()
    iq=np.asarray(Mq@np.ones(q_nodes)).ravel()
    iz=np.asarray(Mz@np.ones(z_nodes)).ravel()
    volume=np.kron(iq,iz)
    end=np.zeros(z_nodes)
    end[[0,-1]]=1
    load=np.kron(iq,end)
    matrix=bmat([[A,csc_matrix(volume[:,None])],
                 [csc_matrix(volume[None,:]),csc_matrix((1,1))]],format="csc")
    solution=spsolve(matrix,np.r_[load,0.])
    u=solution[:-1]
    lam=float(solution[-1])
    residual=float(np.linalg.norm(A@u+volume*lam-load)/np.linalg.norm(load))
    return dict(Q_mm2=-.5*float(load@u),
                volume_residual_mm3=float(volume@u),
                stationarity_relative_residual=residual,
                multiplier_per_mm=lam,
                q_nodes=q_nodes,z_nodes=z_nodes)


def report() -> dict:
    geometry=Channel()
    rows=[]
    for phase in (-1,1):
        state=critical_state(phase,geometry)
        exact=susceptibility_exact(state,geometry.depth_mm,even_modes=320)
        estimates=[susceptibility_fem(state,geometry.depth_mm,n,((n+1)//2))
                   for n in (17,33,65,129)]
        for result in estimates:
            result["absolute_error_mm2"]=abs(result["Q_mm2"]-exact["Q_mm2"])
        rows.append(dict(phase=phase,s_critical_mm=state["s_mm"],
                         kappa_per_mm=state["kappa"],
                         half_arc_mm=state["L"],
                         physical_beta_per_mm=state["beta_wall"],
                         continuum=exact,fem=estimates))
    return dict(status="neutral-plate single-interface cylindrical Jacobi model only",
                theta_side_assumed_deg=geometry.theta_deg,depth_mm=geometry.depth_mm,
                energy_normalization="surface_tension_gamma_times_mm_squared",
                phases=rows,
                difference_Q_mm2=rows[1]["continuum"]["Q_mm2"]-rows[0]["continuum"]["Q_mm2"],
                established=False,experimentally_validated=False,
                limitations=[
                  "Q is not a 3D two-front quartic or a finite-angle nonlinear pressure curve",
                  "All wetting angles and the phase-engineered walls are model inputs, not measured",
                  "No global stability or positive-angle continuation is established",
                  "Finite-mode sum numerically truncated; independently checked by increasing count",
                ])


def csv_content(data: dict) -> str:
    buf=io.StringIO()
    writer=csv.writer(buf,lineterminator="\n")
    writer.writerow(["phase","q_nodes","z_nodes","Q_fem_mm2",
                     "Q_exact_mm2","absolute_error_mm2","stationarity_relative_residual",
                     "volume_residual_mm3"])
    for phase in data["phases"]:
        for r in phase["fem"]:
            writer.writerow([phase["phase"],r["q_nodes"],r["z_nodes"],
                             f'{r["Q_mm2"]:.14g}',f'{phase["continuum"]["Q_mm2"]:.14g}',
                             f'{r["absolute_error_mm2"]:.14g}',
                             f'{r["stationarity_relative_residual"]:.7g}',
                             f'{r["volume_residual_mm3"]:.7g}'])
    return buf.getvalue()


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--csv",action="store_true")
    args=parser.parse_args()
    r=report()
    print(csv_content(r) if args.csv else json.dumps(r,indent=2))
