"""Round 10 extension: independent energy-barrier scaling near the second fold.

Uses the tested full 4-variable, free-contact-angle planar model. Locates
both the local minimum and index-one saddle from equal pressure + volume,
then checks the saddle-node 3/2 barrier law.

Results are MODEL ONLY, not measured jumps or thermally activated rates.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.optimize import brentq
from round10_full_equilibria import PhasePlug

OUT=Path(__file__).parent/"outputs"
OUT.mkdir(exist_ok=True)
SURFACE_TENSION=.072 # illustrative N/m
DEPTH=.38 # mm

def main():
    m=PhasePlug(1,36.5,higher_harmonic_mm=.0004)
    s0,VC=m.critical()
    VF,fold=m.secondary_fold()
    cf=fold["c_um"]/1000
    sf=fold["s_mm"]
    assert abs(VF-.3133490184755)<2e-7
    assert abs(cf-.023600694)<2e-6
    rows=[]
    # A fold has a +stable and -unstable root at identical liquid volume.
    for nL in (.02,.04,.08,.15,.25,.4,.6,1.0,1.5,2.0):
        V=VF+nL/1000
        amp=.0005265*np.sqrt(nL/.005)
        stable=m.station(V,(cf-amp,sf))
        saddle=m.station(V,(cf+amp,sf))
        assert stable["c_um"]<fold["c_um"]<saddle["c_um"]
        assert all(x>0 for x in stable["eigenvalues"])
        assert sum(x < -1e-5 for x in saddle["eigenvalues"])==1
        assert m.caps_clear_walls(stable) and m.caps_clear_walls(saddle)
        dE=saddle["energy_mm"]-stable["energy_mm"]
        assert dE>0
        rows.append(dict(volume_above_fold_nL=nL,
            stable_centre_um=stable["c_um"],
            saddle_centre_um=saddle["c_um"],
            barrier_energy_model_mm=dE,
            illustrative_barrier_pJ=dE*SURFACE_TENSION*DEPTH*1e6,
            scaled_barrier_model_mm_per_nL_3over2=dE/nL**1.5))
    fit=[r for r in rows if r["volume_above_fold_nL"]<=.15]
    x=np.log([r["volume_above_fold_nL"] for r in fit])
    y=np.log([r["barrier_energy_model_mm"] for r in fit])
    exponent=float(np.polyfit(x,y,1)[0])
    assert 1.48<exponent<1.52,exponent
    # The second fold is not the initial symmetry-breaking threshold.
    assert 2.4<1000*(VC-VF)<2.8
    with (OUT/"round10_secondary_fold_barrier.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    summary=dict(scope="planar full Young-angle equilibrium, no empirical data",
       symmetric_onset_volume_uL=VC,
       secondary_fold_volume_uL=VF,
       continuous_asymmetry_interval_nL=1000*(VC-VF),
       secondary_fold_centre_um=fold["c_um"],
       secondary_fold_constrained_eigenvalues=fold["eigenvalues"],
       barrier_exponent_near_fold=exponent,
       notes=["Generic saddle-node exponent, not a novel universal exponent",
              "Not a global minimum-energy-path proof",
              "Stationary Young-angle, equal-pressure, volume and full Hessian checked",
              "No measured wall trace or contact-line trajectories"],
       points=rows)
    (OUT/"round10_secondary_fold_barrier.json").write_text(json.dumps(summary,indent=2))
    print("PASS genuine free-contact-line fold, index-one saddle, barrier 3/2 law")
    print("continuous volume interval nL",summary["continuous_asymmetry_interval_nL"])
    print("near-fold exponent",exponent)


if __name__=="__main__":
    main()
