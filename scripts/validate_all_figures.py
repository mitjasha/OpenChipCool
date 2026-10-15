#!/usr/bin/env python3
import argparse, json
from pathlib import Path
from openchipcool.cfd.validation import validate_export
p=argparse.ArgumentParser()
p.add_argument('--cfd-export',type=Path,required=True)
p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
p.add_argument('--output',type=Path,default=Path('results/cfd/figures4_5_7_report.json'))
a=p.parse_args()
report=validate_export(a.cfd_export,a.root)
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print('Evidence comparison only, not an experimental validation certificate:',a.output)
