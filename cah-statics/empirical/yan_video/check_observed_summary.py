"""Small independent arithmetic and dataset integrity audit for Round 8."""
import csv,math
from pathlib import Path
f=Path(__file__).parent/"observed_video_summary.csv"
rows=list(csv.DictReader(f.open()))
assert len(rows)==2
Q=.01*1000/3600
V=800*800*120/1e9
for r in rows:
    t=float(r["duration_s"])
    assert 0<t<30
    assert math.isclose(float(r["Q1_uL_per_s"]),Q,rel_tol=1e-12)
    assert math.isclose(float(r["nominal_pore_uL"]),V,rel_tol=1e-12)
    assert math.isclose(float(r["pump_during_interval_uL"]),Q*t,rel_tol=1e-11)
    assert math.isclose(float(r["nominal_fill_rate_over_pump"]),.8*V/(Q*t),rel_tol=1e-10)
    assert float(r["duration_mask_min_s"])<=t<=float(r["duration_mask_max_s"])
t60,t120=[float(r["duration_s"]) for r in rows]
print(f"PASS 2 independent CC-BY image stacks; 10-90%% ratio: {t60/t120:.3f} times faster for 120-degree case")
print("Measurements are relative optical signals; NOT independently calibrated liquid volume or burst-pressure predictions.")
