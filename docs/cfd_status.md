# CFD evidence status — OpenChipCool v0.3

| Requirement | Evidence present | Status |
|---|---|---|
| Qu & Mudawar 2002 primary paper | `docs/qu_mudawar_openfoam_protocol.md` | SOURCE VERIFIED |
| 21 full-length microchannels, nominal geometry | `src/openchipcool/cfd/geometry.py` | GENERATED |
| Block-face topology/volume bookkeeping | `mesh_audit.json`, `tests/test_cfd_case.py` | PYTHON PASS |
| OpenFOAM `blockMesh` and `checkMesh` | OpenFOAM not installed in this runtime | NOT RUN |
| OpenFOAM actual CHT solution | No solver in runtime | NOT RUN |
| Experimental pressure comparison | Fig. 4 hand-read markers, no CFD points | NOT VALIDATED |
| Four thermocouple temperatures | Fig. 7 sources marked, not digitized | NOT VALIDATED |
| 3-mesh solved CFD GCI | 3 mesh resolutions specified, none solved | NOT RUN |
| 36-mm new cold plate | 108 channel core prepared | GENERATED ONLY |
| CFD HTC and Δp into FEM/RC/Modelica | No valid CHT data available | NOT INTEGRATED |

Plots in `results/cfd` are **literature/analytical only**, not solver comparisons. Never silently promote `GENERATED` to `VALIDATED`.
