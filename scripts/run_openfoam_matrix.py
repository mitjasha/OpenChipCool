#!/usr/bin/env python3
"""Explicit execution gate for the OpenFOAM.com matrix: never fabricate success.

Use on an OpenFOAM.com installation after sourcing OpenFOAM bashrc.
One case at a time, logs and mesh checks preserved. GCI requires fine/medium/
coarse solver outputs and cannot be computed by this runner alone.
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path
from openchipcool.cfd.foam_setup import write_fields_after_split
from openchipcool.cfd.geometry import geometry_profile
from openchipcool.cfd.reference import experimental_condition

REQUIRED=('blockMesh','splitMeshRegions','checkMesh','chtMultiRegionSimpleFoam')

def shell(cmd,cwd,log):
    with Path(log).open('w') as f:
        f.write('COMMAND: '+' '.join(cmd)+'\n')
        f.flush()
        result=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{' '.join(cmd)} exited {result.returncode}; see {log}")
    return Path(log).read_text(errors='replace')

def run_case(case:Path):
    case=case.resolve()
    info=json.loads((case/'setup_status.json').read_text())
    if info['stage']!='CASE_GENERATED_NOT_SOLVED':
        raise RuntimeError('Refusing reuse of already touched case: generate clean case directory')
    for tool in REQUIRED:
        if not shutil.which(tool):
            raise RuntimeError(f'{tool} unavailable; cannot verify or claim simulation')
    info['stage']='EXECUTING_UNVERIFIED'
    (case/'setup_status.json').write_text(json.dumps(info,indent=2)+'\n')
    try:
        shell(['blockMesh'],case,case/'log.blockMesh')
        shell(['splitMeshRegions','-cellZonesOnly','-overwrite'],case,case/'log.split')
        for r in ('water','copper'):
            out=shell(['checkMesh','-region',r,'-allTopology','-allGeometry'],case,case/f'log.checkMesh.{r}')
            if 'Mesh OK' not in out:
                raise RuntimeError(f'checkMesh {r}: did not report Mesh OK')
        geom=geometry_profile('paper',info['level'])
        c=info['condition']
        write_fields_after_split(case,c['bulk_velocity_m_s'],geom.tin_K,c['qeff_W_cm2'])
        solver_log=shell(['chtMultiRegionSimpleFoam'],case,case/'log.cht')
        if '\nEnd\n' not in solver_log and not solver_log.rstrip().endswith('End'):
            raise RuntimeError('Solver did not print End')
        info['stage']='SOLVER_EXITED_SUCCESSFULLY_NEEDS_CONVERGENCE_AND_VALIDATION'
        info['cfdsimulation_completed']=True
        info['experimentally_validated']=False
    except Exception as ex:
        info['stage']='FAILED_NO_VALID_CFD_EVIDENCE'
        info['exception']=str(ex)
        raise
    finally:
        (case/'setup_status.json').write_text(json.dumps(info,indent=2)+'\n')
    return info

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('case',type=Path)
    ap.add_argument('--execute',action='store_true',help='required to actually run commands')
    args=ap.parse_args()
    if not args.execute:
        print('DRY RUN: explicit --execute required; nothing was solved')
        sys.exit(2)
    try:
        result=run_case(args.case)
        print(result['stage'])
    except Exception as ex:
        print('CFD NOT VALIDATED:',ex,file=sys.stderr)
        sys.exit(2)
