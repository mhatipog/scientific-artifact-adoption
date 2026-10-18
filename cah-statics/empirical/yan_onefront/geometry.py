"""Independent, deliberately restricted corner-rotation calculation.

Never interpret this planar contribution as a measured full 3D burst pressure.
Inputs: reported rounded corner radius ~5.6 micrometres; nominal theta values;
150 micrometre assumed throat width; water-air gamma illustrative 0.072 N/m.
"""
from pathlib import Path
from math import radians,sin,cos,pi
import csv
OUT=Path(__file__).resolve().parent/'outputs'
OUT.mkdir(exist_ok=True)
GAMMA=.072; R=5.6e-6; W=150e-6
rows=[]
for theta in [60,95,120]:
    phi=90-theta if theta<90 else None
    arc=R*radians(phi)*1e6 if phi is not None else None
    dx=R*sin(radians(phi))*1e6 if phi is not None else None
    p0=2*GAMMA*cos(radians(theta))/W/1000
    p90=2*GAMMA*cos(radians(theta+90))/W/1000
    rows.append(dict(nominal_angle_deg=theta,reported_corner_radius_um=5.6,
        assumed_local_width_um=150,
        pressure_start_component_kPa=round(p0,6),
        pressure_after_90deg_turn_component_kPa=round(p90,6),
        predicted_zero_curvature_wall_turn_deg=phi,
        predicted_corner_arclength_um=round(arc,6) if arc is not None else None,
        predicted_projected_displacement_um=round(dx,6) if dx is not None else None))
assert abs(rows[0]['predicted_corner_arclength_um']-2.932153)<1e-5
assert rows[0]['predicted_projected_displacement_um']==2.8
assert rows[1]['predicted_zero_curvature_wall_turn_deg'] is None
with (OUT/'yan_corner_geometric_control.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('PASS geometrical corner sign reversal and dimensional checks')
print('NOT validated: burst pressure or synced event position; nominal angle != advancing angle.')
