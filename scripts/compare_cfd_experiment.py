#!/usr/bin/env python3
"""Fail closed without genuine independent OpenFOAM CFD exported results."""
import argparse, json
from pathlib import Path
from openchipcool.cfd.reference import compare_against_paper

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--cfd-results',type=Path,required=True)
 p.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
 p.add_argument('--out',type=Path,default=Path('results/cfd/comparison.json'))
 a=p.parse_args()
 rows=compare_against_paper(a.cfd_results,a.repo_root)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(dict(status='comparison_completed_NOT_validation_certified',points=rows),indent=2)+'\n')
 print('Wrote',a.out)
 print('Data: Fig. 4 hand-digitized with uncertainties; verify provenance independently')
if __name__=='__main__': main()
