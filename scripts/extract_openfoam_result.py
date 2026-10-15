#!/usr/bin/env python3
"""Convert OpenFOAM v2312 CHT function-object evidence to comparison CSV.

No interpolation of absent data; for the experimental Δp tap-to-tap value,
a separately documented plenum pressure-loss correction is mandatory.
"""
import argparse, csv, json, math
from pathlib import Path
from openchipcool.cfd.reference import experimental_condition
from openchipcool.cfd.geometry import geometry_profile

def last_vector_from_foam(path:Path)->list[float]:
    """Read final numeric line from a probe or surfaceFieldValue ASCII file."""
    data=[]
    for line in path.read_text(errors='replace').splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts=line.replace('(',' ').replace(')',' ').split()
        try:
            row=[float(v) for v in parts]
        except ValueError: continue
        if row and all(math.isfinite(x) for x in row): data.append(row)
    if not data: raise ValueError(f'No numeric time samples: {path}')
    return data[-1]

def latest(root:Path,obj:str,field:str):
    p=root/'postProcessing'/obj
    if not p.is_dir(): raise FileNotFoundError(f'Missing CHT postprocessing object {p}')
    paths=list(p.glob(f'*/{field}'))
    if not paths: raise FileNotFoundError(f'Missing CHT postprocessing field {field} for {obj} under {p}')
    return max(paths,key=lambda f:float(f.parent.name))

def extract(case:Path,plenum_correction:float,plenum_justification:Path):
    case=case.resolve()
    status=json.loads((case/'setup_status.json').read_text())
    if status.get('stage')!='SOLVER_EXITED_SUCCESSFULLY_NEEDS_CONVERGENCE_AND_VALIDATION':
        raise ValueError('OpenFOAM run has not produced a successful solver exit')
    if plenum_correction < 0 or not math.isfinite(plenum_correction):
        raise ValueError('Plenum correction cannot be negative')
    if not plenum_justification.is_file() or len(plenum_justification.read_text().strip())<80:
        raise ValueError('Independent plenum correction justification missing/incomplete')
    values={}
    names={'p_in':'inlet_p_areaAverage','p_out':'outlet_p_areaAverage',
           'mass_in':'inlet_mass_flow','mass_out':'outlet_mass_flow',
           'T_out':'outlet_T_mass_weighted'}
    for tag,obj in names.items():
        last=last_vector_from_foam(latest(case,obj,'surfaceFieldValue.dat'))
        if len(last)<2: raise ValueError('Malformed '+obj)
        values[tag]=last[1]
    tc=last_vector_from_foam(latest(case,'tc_copper','T'))
    if len(tc)!=5: raise ValueError('Expected time and four thermocouple temperatures')
    tin=status['condition']['Tin_C']; q=status['condition']['qeff_W_cm2']
    meta=status['condition']; m_dot=meta['massflow_kg_s']
    mdot_err=abs(abs(values['mass_out'])-abs(values['mass_in']))/max(abs(values['mass_in']),1e-15)
    # Gas/liquid CHT solver phis are mass fluxes. Compare with imposed mdot too.
    mass_in_err=abs(abs(values['mass_in'])-m_dot)/m_dot
    energy=(abs(values['mass_out'])*4180*(values['T_out']-273.15-tin))
    energy_err=abs(energy-meta['power_W'])/meta['power_W']
    if mdot_err>=.001 or mass_in_err>=.01 or energy_err>=.01:
        raise ValueError(f'Mass/energy closure not verified: mnet={mdot_err}, min={mass_in_err}, energy={energy_err}')
    result={'case_id':case.name,'mesh':status['level'], 'Re':meta['Re'],
            'qeff_W_cm2':q,'dp_core_Pa':values['p_in']-values['p_out'],
            'dp_plenums_Pa':plenum_correction,'Tin_C':tin,
            'Tout_massweighted_C':values['T_out']-273.15,
            'mass_error_frac':mdot_err,'energy_error_frac':energy_err,
            'solver_log_path':str(case/'log.cht'),
            'mesh_water_check_path':str(case/'log.checkMesh.water'),
            'mesh_copper_check_path':str(case/'log.checkMesh.copper'),
            'plenum_justification_file':str(plenum_justification.resolve()),
            **{f'Ttc{i}_rise_K':x-273.15-tin for i,x in enumerate(tc[1:],1)}}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('case',type=Path)
    p.add_argument('--plenum-correction-pa',type=float,required=True)
    p.add_argument('--plenum-justification',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    vals=extract(a.case,a.plenum_correction_pa,a.plenum_justification)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=vals.keys());w.writeheader();w.writerow(vals)
    print('Evidence CSV written. NOT EXPERIMENTALLY VALIDATED.',a.out)
