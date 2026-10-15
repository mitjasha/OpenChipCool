"""Hydraulic analytical cross-check and honest, figure-digitized evidence comparison.

No function in this module runs or impersonates OpenFOAM.  All experimental
pressure markers except explicitly documented textual anchor are approximate
hand readings from Qu & Mudawar (2002) Fig. 4(a,b).
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import csv
import math
from .geometry import Geometry, geometry_profile, geometry_metrics

DOI='10.1016/S0017-9310(01)00337-4'
PAPER='Qu & Mudawar (2002), IJHMT 45, 2549-2565'

@dataclass(frozen=True)
class Water:
    rho_kg_m3: float=999.1 # 15 deg C, approximate constant property
    mu_Pa_s: float=1.138e-3
    cp_J_kgK: float=4180.
    k_W_mK: float=0.589


def experimental_condition(Re:float, qeff_W_cm2:float=100., g:Geometry|None=None, water:Water|None=None)->dict:
    g=g or geometry_profile()
    w=water or (Water(rho_kg_m3=995.7,mu_Pa_s=0.000798,cp_J_kgK=4178.,k_W_mK=.613)
                if g.name=='openchipcool_36mm_concept' else Water())
    g.validate()
    if Re<=0 or qeff_W_cm2<=0:
        raise ValueError('Re and nominal imposed heat flux must be >0')
    dh=geometry_metrics(g)['hydraulic_diameter_m']
    u=Re*w.mu_Pa_s/(w.rho_kg_m3*dh)
    V=u*g.n_channels*g.channel_width_m*g.channel_depth_m
    Q=qeff_W_cm2*1e4*g.length_m*g.width_m
    mdot=w.rho_kg_m3*V
    return dict(Re=float(Re),qeff_W_cm2=qeff_W_cm2,
                bulk_velocity_m_s=u,flow_L_min=V*60000.,
                massflow_kg_s=mdot,power_W=Q,
                adiabatic_outlet_rise_K=Q/(mdot*w.cp_J_kgK),
                Tin_C=g.tin_K-273.15,
                outlet_energy_upper_C=g.tin_K-273.15+Q/(mdot*w.cp_J_kgK))


def laminar_core_pressure_pa(Re:float,g:Geometry|None=None,w:Water|None=None)->float:
    """Hagen-Poiseuille f_D Re approximation, developed *channel only*.

    Shah/London rectangular-channel laminar f_D*Re polynomial. This is not
    the measured plenum-to-plenum pressure drop of the paper (Eq. 18).
    """
    g=g or geometry_profile()
    w=w or (Water(rho_kg_m3=995.7,mu_Pa_s=.000798,cp_J_kgK=4178.,k_W_mK=.613)
            if g.name=='openchipcool_36mm_concept' else Water())
    c=experimental_condition(Re,g=g,water=w)
    a=min(g.channel_width_m,g.channel_depth_m)/max(g.channel_width_m,g.channel_depth_m)
    po=96*(1-1.3553*a+1.9467*a*a-1.7012*a**3+0.9564*a**4-0.2537*a**5)
    return (po/Re)*(g.length_m/geometry_metrics(g)['hydraulic_diameter_m'])*0.5*w.rho_kg_m3*c['bulk_velocity_m_s']**2


def _digitized(root:Path, flux:float)->list[dict]:
    path=Path(root)/'data/literature/qu_mudawar_2002_pressure_approx.csv'
    with path.open(newline='') as f: data=list(csv.DictReader(f))
    return [dict(Re=float(d['Re']),dp_bar=float(d['dp_bar_visual_approx']))
            for d in data if float(d.get('qeff_W_cm2',100))==float(flux)]


def compare_against_paper(prediction_csv:Path, repo_root:Path):
    """Compare actual exported CFD CSV against approximate *published* markers.

    Nothing is treated as a measurement of this project's own cold plate. To
    prevent accidental false reporting, require an independent CFD CSV file;
    measured/analytical curves are NOT silently substituted for CFD output.
    """
    p=Path(prediction_csv)
    if not p.is_file():
        raise FileNotFoundError('No CFD export found: OpenFOAM VALIDATION NOT RUN')
    with p.open(newline='') as f:
        rows=list(csv.DictReader(f))
    needed={'Re','qeff_W_cm2','dp_core_Pa','dp_total_Pa','Tout_C','energy_error_fraction'}
    if not rows or not needed.issubset(rows[0]):
        raise ValueError('CFD export must provide '+', '.join(sorted(needed)))
    results=[]
    for r in rows:
        Re=float(r['Re']); flux=float(r['qeff_W_cm2']); dp=float(r['dp_total_Pa'])
        core=float(r['dp_core_Pa']); eps=float(r['energy_error_fraction'])
        samples=_digitized(repo_root,flux)
        if len(samples)<2 or not (samples[0]['Re']<=Re<=samples[-1]['Re']):
            raise ValueError('Extrapolation outside source-marker Re range is not an experimental comparison')
        import numpy as np
        exp_pa=float(np.interp(Re,[s['Re'] for s in samples],[s['dp_bar'] for s in samples]))*1e5
        results.append(dict(Re=Re,qeff_W_cm2=flux,CFD_core_Pa=core,
                           CFD_plenum_corrected_Pa=dp,exp_fig4_approx_Pa=exp_pa,
                           pressure_error_fraction=(dp-exp_pa)/exp_pa,
                           Tout_C=float(r['Tout_C']),energy_error_fraction=eps,
                           comparison_type='paper Fig4 visually approximated; not raw experimental data'))
    return results
