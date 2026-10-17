"""Numerical pressure-volume continuation for full 3D free-boundary capillary surfaces.

No imposed meniscus pressure, no fitted experimental jump data. The surfaces
are reflection-symmetric in y and z; stability is local to that class.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.polynomial import polyfit, polyder, polyval
from surface_3d import solve

OUT=Path(__file__).parent
ang_grid=[90,85,80]
scrit={-1:0.7003307,1:0.697477}
o=[]
for sign in [-1,1]:
    b=sign*.0005
    sc=scrit[sign]
    for plate in ang_grid:
        for j,s in enumerate(np.linspace(sc-.012,sc+.012,19)):
            r=solve(float(s),plate_deg=plate,b=b,deg_e=4,deg_z=3,ne=32,nz=10)
            assert r['opt'].success and r['error']<1.5e-5,(s,r['error'],r['opt'].message)
            o.append(dict(phase=sign,plate_angle_deg=plate,sref_mm=s,pressure_invmm=r['p'],single_front_volume_uL=r['volume'],
                          stationarity_normalized=r['error'],correction_max_um=1000*max(abs(r['c']))))
        print('FINISHED',sign,plate,'count',len(o),flush=True)
with (OUT/'pressure_volume_3d.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(o[0].keys()));wr.writeheader();wr.writerows(o)

summ=[]
for sign in [-1,1]:
 for plate in ang_grid:
  dat=[r for r in o if r['phase']==sign and r['plate_angle_deg']==plate]
  dat=sorted(dat,key=lambda r:r['single_front_volume_uL'])
  q=np.array([r['single_front_volume_uL'] for r in dat]); p=np.array([r['pressure_invmm'] for r in dat]);
  q0=q[9];scale=(q.max()-q.min())/2
  xx=(q-q0)/scale
  for order in [4,5,6]:
   coeff=polyfit(xx,p,order)
   d1=polyder(coeff);d2=polyder(d1);d3=polyder(d2)
   roots=np.roots(np.array(d1[::-1]))
   real=[float(x.real) for x in roots if abs(x.imag)<1e-5 and -.8<x.real<.8]
   if not real:print('no extremum',sign,plate,order);continue
   root=min(real,key=abs)
   p2=polyval(root,d2)/scale**2
   p3=polyval(root,d3)/scale**3
   pred='supercritical' if p3<0 else 'subcritical'
   rms=np.sqrt(np.mean((polyval(xx,coeff)-p)**2))
   summ.append(dict(phase=sign,plate_angle_deg=plate,fit_degree=order,critical_front_volume_uL=q0+root*scale,
      critical_pressure_invmm=polyval(root,coeff),p2_invmm_per_uL2=p2,p3_invmm_per_uL3=p3,
      local_onset_type=pred,pressure_fit_rms_invmm=rms,span_uL=2*scale))
   print('FIT',sign,plate,order,'threshold',q0+root*scale,'p3',p3,'class',pred,'RMSE',rms)
with (OUT/'normal_form_3d.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(summ[0].keys()));wr.writeheader();wr.writerows(summ)