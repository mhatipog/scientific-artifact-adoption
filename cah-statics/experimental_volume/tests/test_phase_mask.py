"""Network-free regression tests. No original TIFF or Actions runtime needed."""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from segment_original_tiffs import find_pore, crossing, PORE_VOLUME_UL, PUMP_UL_S


def main():
    assert abs(PORE_VOLUME_UL - .0768) < 1e-12
    assert abs(PUMP_UL_S - .002777777777777778) < 1e-12
    assert abs(crossing(np.array([0.,1.,2.]), np.array([0.,.5,1.]),.75)-1.5)<1e-12
    assert crossing(np.array([0.,1.]), np.array([0.,.5]),.9) is None

    # Both observed source offsets independently reconstructed on synthetic
    # dark-wall channels; this is a deterministic geometry test, not raw data.
    for yy in (31,62):
        a=np.full((512,512),3300,dtype=np.float32)
        x0,x1=96,418
        a[yy-1:yy+2,x0:x1]=1000
        a[yy+321:yy+324,x0:x1]=1000
        a[yy:yy+323,x0-1:x0+2]=1000
        a[yy:yy+323,x1-1:x1+2]=1000
        assert find_pore(a)==(96,yy,418,yy+322)

    summary=json.loads((HERE/"derived"/"measurements.json").read_text())
    with (HERE/"derived"/"segmented_saturation_10hz.csv").open(newline="") as f:
        rows=list(csv.DictReader(f))
    for deg,min_n,terminal in ((60,125,.93967),(120,215,.91232)):
        r=[v for v in rows if int(v["angle_deg"])==deg]
        assert len(r)>=min_n
        assert float(r[0]["clip_t_s"])==0.
        assert all(float(v["mask_lo"])>=0 for v in r)
        assert abs(float(r[-1]["mask_reference"])-terminal)<.003
        output=summary["reference_mask_calculation"][str(deg)]
        assert .045<output["net_storage_uL"]<.056
        assert 4<output["flow_ratio"]<16
    print("PASS network-free geometry, crossing, source completeness and volume-budget regressions")


if __name__ == "__main__":
    main()
