"""Published-figure reading ledger; no raw camera/pressure data downloaded.

Reproduces only geometric arithmetic, not Yan et al. pressure traces.
"""
from math import pi
import csv
from pathlib import Path

pore_uL = 800*800*120/1e9
pump_uL_per_s = .01*1000/3600
print("Pore volume [uL]:", pore_uL)
print("Pump-only filling time [s]:", pore_uL/pump_uL_per_s)
assert abs(pore_uL-.0768) < 1e-12
assert abs(pore_uL/pump_uL_per_s-27.648) < 1e-8
cases=[]
for theta, duration, s0,s1,dp in [
    (60,5.5,.10,1.00,300), (95,.5,.10,.95,430), (120,2.5,.10,.92,450)]:
    delta=(s1-s0)*pore_uL
    displaced=duration*pump_uL_per_s
    cases.append(dict(angle_deg=theta, rough_dt_s=duration,
                      approximate_fill_ratio=delta/displaced,
                      inferred_C_uL_per_kPa=(delta-displaced)/dp*1000))
out=Path(__file__).parent/"outputs"
out.mkdir(exist_ok=True)
with (out/"rough_mass_balance.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=cases[0].keys())
    writer.writeheader()
    writer.writerows(cases)
print("CAUTION: estimated from separate published figures, NOT raw data.",cases)
