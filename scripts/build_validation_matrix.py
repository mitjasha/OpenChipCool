#!/usr/bin/env python3
"""Create source-traceable three-grid CHT input matrix without claiming runs."""
import argparse, json, sys
from pathlib import Path
from openchipcool.cfd.geometry import geometry_profile,write_geometry_case
from openchipcool.cfd.foam_setup import write_fv_dictionaries,write_thermophysics
from openchipcool.cfd.reference import experimental_condition

# Fig 7 target Reynolds numbers, 100 and 200 W/cm2, three grids each.
OPERATING_POINTS=((100.,226.),(100.,890.),(100.,1454.),
                  (200.,442.),(200.,864.),(200.,1220.))

def build(out:Path,levels=('coarse','medium','fine')):
    records=[]
    for q,re in OPERATING_POINTS:
        for level in levels:
            g=geometry_profile('paper',level)
            name=f'qm_q{int(q)}_Re{int(re)}_{level}'
            case=out/name
            metrics=write_geometry_case(g,case)
            ex=write_fv_dictionaries(case,g,re,q)
            write_thermophysics(case,g)
            status={'stage':'CASE_GENERATED_NOT_SOLVED','cfdsimulation_completed':False,
                    'experimentally_validated':False,'condition':ex,
                    'planned_hex_cells':metrics['cells_total'],'level':level,
                    'plenum_geometry':'NOT PRESENT','cover_conduction':'NOT PRESENT',
                    'temperature_dependent_viscosity':'NOT IMPLEMENTED',
                    'tc_holes_and_slots':'NOT PRESENT'}
            (case/'setup_status.json').write_text(json.dumps(status,indent=2)+'\n')
            records.append({'case_id':name,'mesh':level,'Re':re,'qeff_W_cm2':q,
                            'cells_planned':metrics['cells_total'],
                            'planned_flow_L_min':ex['flow_L_min'],
                            'power_W':ex['power_W'],'cfd_status':status['stage']})
    manifest={'paper':'https://doi.org/10.1016/S0017-9310(01)00337-4',
              'count':len(records),'source_figure_numbers':[4,5,7],
              'status':'INPUTS_ONLY_NO_SIMULATION',
              'cases':records}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=Path('cfd_cases/validation_matrix'))
    ap.add_argument('--levels',nargs='+',choices=('coarse','medium','fine'),default=('coarse','medium','fine'))
    args=ap.parse_args()
    v=build(args.output,args.levels)
    print('Generated',v['count'],'CFD input cases. NO SOLVER RESULTS.',file=sys.stdout)
