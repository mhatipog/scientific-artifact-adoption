import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from numpy.polynomial.polynomial import polyfit,polyder,polyval
from surface_3d import solve,D,K

here=Path(__file__).parent
read=lambda f:list(csv.DictReader(open(here/f)))
data=read('normal_form_3d.csv')
ext=read('normal_form_3d_extended.csv')
conv=read('convergence_onset.csv')+read('convergence_onset_deep.csv')
ref=[[r for r in data if int(r['fit_degree'])==6 and int(r['phase'])==p] for p in [-1,1]]
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,
'axes.titlesize':13,'axes.labelsize':10,'legend.fontsize':9,'figure.facecolor':'#fdfdfd'})
colors={-1:'#c55f4f',1:'#167c8d'}

fig=plt.figure(figsize=(13.7,9.8),constrained_layout=True)
grid=GridSpec(2,2,figure=fig,height_ratios=[1,1],hspace=.25,wspace=.15)
ax=fig.add_subplot(grid[0,0])
for p in [-1,1]:
 lines=ref[0] if p==-1 else ref[1]
 ax.plot([float(r['plate_angle_deg']) for r in lines], [float(r['p3_invmm_per_uL3']) for r in lines],
         '-o',lw=2.8,color=colors[p],label=('A · minus phase' if p==-1 else 'B · plus phase'))
ax.axhline(0,color='#444',lw=1)
ax.set(xlim=(79,91),xlabel='Top and bottom Young angle (degrees)',
 ylabel=r'$d^3p/dQ^3$ (mm$^{-1}$/µL$^3$)',
 title='A | Opposite onset classes survive non-neutral plates')
ax.text(80.3,2800,'Subcritical region: local finite jump possible',color='#a54b3d',fontsize=9)
ax.text(80.3,-2350,'Supercritical region: local smooth asymmetry',color='#167c8d',fontsize=9)
ax.legend(frameon=False,loc='upper left')
ax.grid(alpha=.16)

ax=fig.add_subplot(grid[0,1])
for p in [-1,1]:
 lines=ref[0] if p==-1 else ref[1]
 ax.plot([float(r['plate_angle_deg']) for r in lines], [float(r['critical_pressure_invmm']) for r in lines],
         '-o',lw=2.8,color=colors[p],label=('A · minus phase' if p==-1 else 'B · plus phase'))
ax.set(xlim=(79,91),xlabel='Top and bottom Young angle (degrees)',ylabel='Critical mean curvature p = ΔP/γ (mm⁻¹)',
       title='B | Depth wetting strongly shifts capillary pressure')
ax.grid(alpha=.16)

ax=fig.add_subplot(grid[1,0])
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
for phase in [-1,1]:
 for angle in [90,85,80]:
  row=next(r for r in data if int(r['fit_degree'])==6 and int(r['phase'])==phase and int(r['plate_angle_deg'])==angle)
  ax.scatter(angle,float(row['p3_invmm_per_uL3'])/1000,s=85 if phase==-1 else 70,
             marker='s' if phase==-1 else 'o',color=colors[phase],edgecolors='white',zorder=5)
for row in conv:
 if int(row['fit_degree'])!=6:continue
 deg=row['mesh'].split()[0]
 if int(row['plate_angle_deg'])>70:continue
 x=int(row['plate_angle_deg'])+({'3x2':-.16,'5x4':0.,'6x5':.16}.get(deg,0))
 ax.scatter(x,float(row['p3'])/1000,s=40,marker='x',color='#888',alpha=.8)
ax.axhline(0,lw=.8,color='#777')
ax.axvspan(59,71,alpha=.08,color='#d88b20')
ax.set(xlim=(59,91),xlabel='Plate angle (degrees)',ylabel=r'$p^{(3)}(Q_c)$ (×10³)',
       title='C | Low-angle transition is NOT numerically converged')
ax.text(59.4,2.6,'Grey ×: low-angle results change with resolution',color='#805715',fontsize=9)
ax.grid(alpha=.12)

ax=fig.add_subplot(grid[1,1])
# Draw a concrete meniscus cut across transverse and depth directions, with true 3D numeric solutions.
b=-.0005
r=solve(.7,plate_deg=80,b=b,deg_e=4,deg_z=3,ne=32,nz=10)
from numpy.polynomial.legendre import Legendre
eta=np.linspace(-1,1,240);z=np.linspace(-1,1,100)
R,C,angle=r['ref'];w0=1+.05*np.cos(K*.7)+b*(np.cos(2*K*.7)-np.cos(4*K*.7))
base=C-np.sqrt(R**2-(eta*w0/2)**2)
shape=np.zeros((len(eta),len(z)))
for val,(i,j) in zip(r['c'],r['modes']):shape+=val*np.outer(Legendre.basis(2*i)(eta),Legendre.basis(2*j)(z))
for zz,label in [(0,'mid-depth'),(49,'near top/bottom plate')]:
 ax.plot(eta*w0*500,(base+shape[:,zz])*1000,label=label,lw=2)
ax.plot(eta*w0*500,base*1000,'--',color='#777',label='neutral-plate exact planar arc')
ax.set(xlabel='Transverse coordinate (µm)',ylabel='Axial interface position (µm)',title='D | Three-dimensional interface deformation, θplate=80°')
ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
fig.suptitle('When wetting depth becomes an independent control parameter',fontsize=20,weight='bold')
fig.text(.5,-.012,'Research calculations only · prescribed ideal wall profiles, θside=36.5° · NOT calibrated against experiment',
         ha='center',fontsize=10,color='#555')
fig.savefig(here/'round5_3d_validation.png',dpi=200,bbox_inches='tight')
fig.savefig(here/'round5_3d_validation.svg',bbox_inches='tight')
plt.close(fig)
print('saved',here/'round5_3d_validation.png')