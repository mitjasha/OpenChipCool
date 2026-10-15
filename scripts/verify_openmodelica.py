#!/usr/bin/env python3
"""Actually execute the generated .mo via omc and inspect resulting CSV.

When no omc is installed this exits 2: this is NOT successful validation.
"""
from pathlib import Path
import shutil,subprocess,sys,csv,os
root=Path(__file__).resolve().parents[1]
omc=shutil.which('omc')
if not omc:
    print('UNVERIFIED: omc (OpenModelica) is not installed/in PATH.', file=sys.stderr)
    sys.exit(2)
mo=root/'results'/'OpenChipCoolLoop.mo'
if not mo.exists():
    print('Run: python -m openchipcool.cli benchmark',file=sys.stderr)
    sys.exit(2)
proc=subprocess.run([omc,'modelica/simulate.mos'],cwd=root,capture_output=True,text=True,timeout=180)
print(proc.stdout); print(proc.stderr,file=sys.stderr)
if proc.returncode!=0 or 'Simulation Failed' in proc.stdout or 'Error:' in proc.stdout:
    print('FAIL: OpenModelica compilation/simulation failure',file=sys.stderr)
    sys.exit(1)
res=root/'OpenChipCoolLoop_res.csv'
if not res.exists():
    print('FAIL: expected OpenChipCoolLoop_res.csv missing',file=sys.stderr)
    sys.exit(1)
with res.open(newline='') as f:
    data=list(csv.DictReader(f))
if len(data)<2 or float(data[-1]['time'])<59:
    print('FAIL: simulated trace incomplete',file=sys.stderr);sys.exit(1)
# Cross-check compiled Modelica with Python under same piecewise forcing.
import numpy as np
from openchipcool.rc import simulate_loop,RCParameters
from openchipcool.config import Spec
from pathlib import Path
import json
sys.path.insert(0,str(root/'src'))
from openchipcool.rc import simulate_loop,RCParameters
from openchipcool.config import Spec
spec=Spec()
rc=RCParameters.from_dict(json.loads((root/'results'/'rc_calibration.json').read_text()))
t=np.asarray([float(row['time']) for row in data])
power=np.full((len(t),4),200.)
power[t<5,:]=100.
power[(t>=15)&(t<35),0]=350.
flows=np.where(t<25,4.,.2)
reference=simulate_loop(t,power,flows,rc,spec)
modelica_die=np.array([float(data[-1][f'Tdie[{i}]'])-273.15 for i in range(1,5)])
err=float(np.max(np.abs(modelica_die-reference['die'][-1])))
print(f'Python vs Modelica final die-temperature difference: {err:.3f} K')
if not np.isfinite(err) or err>3.0:
    print('FAIL: Python / OpenModelica end-state mismatch >3 K',file=sys.stderr)
    sys.exit(1)
print('PASS: OpenModelica compiled, simulated, generated',len(data),'timesteps, and agreed with Python ODE within 3 K.')
