#!/usr/bin/env python3
"""Create standalone OpenFOAM multi-region CHT geometry and explicit provenance."""
import argparse, json
from pathlib import Path
from openchipcool.cfd.geometry import geometry_profile,write_geometry_case
from openchipcool.cfd.reference import experimental_condition

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--profile',choices=['paper','design36'],default='paper')
    p.add_argument('--level',choices=['coarse','medium','fine'],default='coarse')
    p.add_argument('--Re',type=float,default=890.)
    p.add_argument('--qeff',type=float,default=100.,help='Heating W/cm2, only for paper mode')
    p.add_argument('--output',type=Path,default=Path('cfd_cases/generated'))
    a=p.parse_args()
    g=geometry_profile(a.profile,a.level)
    meta=write_geometry_case(g,a.output)
    if a.profile=='paper':
        conditions=experimental_condition(a.Re,a.qeff,g)
        (a.output/'experimental_condition.json').write_text(json.dumps(conditions,indent=2)+'\n')
        print(json.dumps(conditions,indent=2))
    print('Wrote:',a.output.resolve())
    print('Channels:',meta['channels'],'blocks:',meta['blocks'],'hex cells:',meta['cells_total'])
    print('STATUS: geometry generated; CFD NOT RUN, comparison NOT VALIDATED')
if __name__=='__main__': main()
