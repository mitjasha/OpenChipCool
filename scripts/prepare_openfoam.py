#!/usr/bin/env python3
"""Generate a CHT case with honest preflight and optional OpenFOAM execution."""
from pathlib import Path
import argparse, json, shutil, subprocess
from openchipcool.cfd.geometry import geometry_profile,write_geometry_case
from openchipcool.cfd.reference import experimental_condition
from openchipcool.cfd.foam_setup import write_fv_dictionaries,write_thermophysics,write_fields_after_split

def checked(cmd,cwd):
    print('RUN:', ' '.join(cmd), flush=True)
    subprocess.run(cmd,cwd=cwd,check=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--profile',choices=['paper','design36'],default='paper')
    p.add_argument('--level',choices=['coarse','medium','fine'],default='coarse')
    p.add_argument('--Re',type=float,default=None,help='Inlet Reynolds number; default 890 for paper')
    p.add_argument('--flow-lpm',type=float,default=None,help='Alternative to --Re; design36 default is 4 L/min')
    p.add_argument('--qeff',type=float,default=None,help='Bottom planform heat flux W/cm²; paper default 100')
    p.add_argument('--power-W',type=float,default=None,help='Total heating power; design36 default 800 W')
    p.add_argument('--case',type=Path,default=Path('cfd_cases/current'))
    p.add_argument('--mesh',action='store_true',help='Require OpenFOAM tools, mesh, split and check regions')
    a=p.parse_args()
    g=geometry_profile(a.profile,a.level)
    if a.Re is not None and a.flow_lpm is not None:
        p.error('Choose --Re OR --flow-lpm')
    if a.qeff is not None and a.power_W is not None:
        p.error('Choose --qeff OR --power-W')
    effective_flux=(a.power_W/(g.length_m*g.width_m)/1e4 if a.power_W is not None
                    else a.qeff if a.qeff is not None
                    else 100. if a.profile=='paper' else 800/(g.length_m*g.width_m)/1e4)
    if a.flow_lpm is not None or (a.profile=='design36' and a.Re is None):
        from openchipcool.cfd.reference import Water
        water=(Water(rho_kg_m3=995.7,mu_Pa_s=.000798,cp_J_kgK=4178.,k_W_mK=.613)
               if a.profile=='design36' else Water())
        flow=(a.flow_lpm if a.flow_lpm is not None else 4.)*1e-3/60
        dh=2*g.channel_width_m*g.channel_depth_m/(g.channel_width_m+g.channel_depth_m)
        u=flow/(g.n_channels*g.channel_width_m*g.channel_depth_m)
        re=u*water.rho_kg_m3*dh/water.mu_Pa_s
    else:
        re=a.Re if a.Re is not None else 890.
    case=a.case.resolve()
    metrics=write_geometry_case(g,case)
    ex=write_fv_dictionaries(case,g,re,effective_flux)
    write_thermophysics(case,g)
    (case/'setup_status.json').write_text(json.dumps(dict(status='geometry_and_dictionary_generation_only',
        NOT_OPENFOAM_VALIDATED=True,mesh_cells_planned=metrics['cells_total'],condition=ex),indent=2)+'\n')
    print(json.dumps({'path':str(case),'planned_cells':metrics['cells_total'],'setup':'GENERATED, NOT SOLVER-VERIFIED'},indent=2))
    if not a.mesh:
        print('Run with --mesh only on a host sourcing OpenFOAM.com (target syntax: v2312, NOT locally solver-tested).')
        return
    tools=('blockMesh','splitMeshRegions','checkMesh','chtMultiRegionSimpleFoam')
    absent=[k for k in tools if not shutil.which(k)]
    if absent:
        raise SystemExit('FAIL: OpenFOAM commands missing: '+', '.join(absent)+'. No CHT claim.')
    checked(['blockMesh'],case)
    checked(['splitMeshRegions','-cellZonesOnly','-overwrite'],case)
    checked(['checkMesh','-region','water','-allTopology','-allGeometry'],case)
    checked(['checkMesh','-region','copper','-allTopology','-allGeometry'],case)
    write_fields_after_split(case,ex['bulk_velocity_m_s'],g.tin_K,effective_flux)
    (case/'setup_status.json').write_text(json.dumps(dict(status='meshed_checked_not_solved',
        NOT_OPENFOAM_VALIDATED=True,geometry=metrics,condition=ex),indent=2)+'\n')
    print('Geometry check passed. To solve (UNVERIFIED SETTINGS): cd',case,'&& chtMultiRegionSimpleFoam > log.cht 2>&1')

if __name__=='__main__':main()
