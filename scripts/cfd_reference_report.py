#!/usr/bin/env python3
"""Generate PAPER-vs-ANALYTICAL screening chart. Intentionally no CFD curve."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib.pyplot as plt
from openchipcool.cfd.geometry import geometry_profile,geometry_metrics
from openchipcool.cfd.reference import experimental_condition,laminar_core_pressure_pa


def main():
 root=Path(__file__).resolve().parents[1]
 out=root/'results/cfd';out.mkdir(parents=True,exist_ok=True)
 g=geometry_profile('paper','coarse')
 rows=list(csv.DictReader((root/'data/literature/qu_mudawar_2002_pressure_approx.csv').open()))
 re=np.linspace(139,1672,150)
 dp=np.array([laminar_core_pressure_pa(v,g)/1e5 for v in re])
 fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
 colors=['tab:blue','tab:orange']
 for i,q in enumerate((100,200)):
  r=np.array([float(d['Re']) for d in rows if float(d['qeff_W_cm2'])==q])
  p=np.array([float(d['dp_bar_visual_approx']) for d in rows if float(d['qeff_W_cm2'])==q])
  ax.scatter(r,p,label=f'Qu & Mudawar Fig. 4(a) ~read, {q} W/cm²',color=colors[i],s=24,alpha=.85)
 ax.plot(re,dp,'k--',lw=2,label='Analytic developed laminar CORE only, $\\mu(15°C)$')
 ax.set(xlabel='Reynolds number $Re_{in}$',ylabel='Pressure drop (bar)',
        title='Screening baseline — NOT OpenFOAM results, NOT experimental validation')
 ax.grid(alpha=.25);ax.legend(fontsize=8)
 fig.savefig(out/'pressure_reference_analytical_only.png',dpi=170);plt.close(fig)
 points=[]
 for R in (139,226,442,890,1220,1672):
  a=experimental_condition(R,g=g)
  points.append(dict(Re=R,flow_L_min=a['flow_L_min'],Q_W=a['power_W'],
                     core_laminar_delta_p_bar=laminar_core_pressure_pa(R,g)/1e5,
                     adiabatic_temp_rise_K=a['adiabatic_outlet_rise_K']))
 with (out/'analytic_baseline.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=points[0].keys());w.writeheader();w.writerows(points)
 (out/'status.json').write_text(json.dumps({
   'status':'GEOMETRY_AND_ANALYTIC_REFERENCE_ONLY',
   'openfoam_executed':False,'cht_executed':False,'validated':False,
   'source':'Qu & Mudawar 2002 Figure 4 visually approximated (except 0.86 bar at Re1672 textual)',
   'warning':'Paper measures plenum-to-plenum pressure; analytical curve covers channel core only',
   'paper_case_mesh_cells':geometry_metrics(g)['cells_total'],
   'new_36mm_plate_mesh_cells':geometry_metrics(geometry_profile('design36','coarse'))['cells_total'],
   'comparison_required':'3 actual solved meshes x 2 heat fluxes x multiple Re; CFD flow and energy and temperatures'
 },indent=2)+'\n')
 print(out/'pressure_reference_analytical_only.png')

if __name__=='__main__':main()
