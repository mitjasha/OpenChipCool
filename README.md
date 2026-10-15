# OpenChipCool

**Multiscale thermal modeling for liquid-cooled chiplet packages.**

OpenChipCool connects a 3D finite-element heat-conduction model, a compact
thermal RC network and a simplified liquid-cooling loop. It explores how
uneven chip power, workload changes and reduced coolant flow affect the
temperatures of a four-chiplet package.

The project includes Python simulation tools, a Streamlit dashboard,
Modelica export and parametric OpenFOAM case generation for microchannel
cold plates. Current version: **0.4.0**.

![Steady 3D FEM temperature distribution](results/temperature_field.png)

## Capabilities

| Component | What it does |
|---|---|
| 3D thermal FEM | Steady and transient conduction through silicon, TIM, epoxy and copper using linear tetrahedral elements |
| Chip power maps | Four independently powered dies with spatially concentrated hotspots |
| Reduced thermal model | Fits an eight-node RC network to FEM transients and evaluates it on a separate power sequence |
| Cooling loop | Couples the RC network to a well-mixed coolant volume and a constant-UA radiator |
| Interactive dashboard | Explores load steps, chip power pulses, power surges and pump-flow reduction; exports CSV traces |
| Modelica export | Generates a standalone lumped loop using the same RC coefficients as Python |
| CFD preparation | Generates copper/water OpenFOAM cases for a 21-channel literature benchmark and a 108-channel, 36 mm plate |
| Numerical checks | Energy balances, geometry-aligned mesh refinement, conditional GCI estimates and FEM/RC hold-out comparison |

![FEM and RC temperatures on a separate power sequence](results/benchmark.png)

## Quick start

Requires **Python 3.10 or newer**. OpenFOAM and OpenModelica are optional;
the Python benchmark and dashboard run independently of them.

```bash
git clone https://github.com/mitjasha/OpenChipCool.git
cd OpenChipCool
python -m venv .venv
```

Activate the environment on Linux or macOS:

```bash
source .venv/bin/activate
```

Or in Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Install the project and run the benchmark and dashboard:

```bash
python -m pip install -e ".[dashboard,dev]"
python -m openchipcool.cli benchmark --refinement 2 --output results
streamlit run dashboard.py
```

The benchmark generates the RC calibration consumed by the dashboard.
Before calibration, the dashboard uses parameters estimated from geometry.
For a smaller initial run, use `--refinement 1`.

Run the Python checks and mesh study separately:

```bash
python -m pytest -q
python -m openchipcool.cli mesh-study --refinements 1,2,3 --output results
```

## Model structure

```text
Four spatial die power maps
            |
3D tetrahedral thermal FEM ---- mesh and energy checks
            |
Four die-average temperature histories
            |
Eight-node thermal RC network ---- separate hold-out comparison
            |
     +------+------+
     |             |
Python coolant   Modelica coolant
loop            loop (.mo)
     |
Streamlit dashboard and interactive HTML report
```

The FEM solver uses NumPy and SciPy. The reduced model contains four die
nodes and four copper-quadrant nodes; coupling it to the coolant adds a
ninth thermal state. The dashboard runs this reduced model, while the 3D
calculation and RC fitting run through the benchmark command.

## Reference configuration

| Parameter | Value |
|---|---:|
| Copper base | 36 × 36 × 3 mm |
| Silicon dies | Four, each 12 × 12 × 0.5 mm; 14 mm pitch |
| Thermal interface | 0.1 mm; conductivity 5 W/(m·K) |
| Epoxy filler conductivity | 0.7 W/(m·K) |
| Silicon conductivity | 130 W/(m·K) |
| Reference power | 800 W total |
| Reference coolant inlet | 30 °C |
| Reference flow | 4 L/min |
| Coolant volume | 10 L equivalent water capacity |
| Radiator conductance | 300 W/K |
| Illustrative temperature threshold | 95 °C |

This is a synthetic package configuration. Its equivalent cold-plate
convection coefficient is derived from approximate readings of published
CFD Nusselt curves in Qu & Mudawar (2002), with idealized fin and coolant
conditions. It is not a measured coefficient for the 36 mm plate.

## Generated outputs

The benchmark writes its outputs to `results/` by default:

| File | Contents |
|---|---|
| `benchmark.png`, `holdout_comparison.csv` | FEM/RC comparison on a power sequence separate from training |
| `temperature_field.png` | Die-top view and a vertical FEM temperature slice |
| `steady_temperature_points.csv`, `steady_temperature.vtk` | Nodal temperature field; the VTK file opens in ParaView |
| `benchmark_summary.json`, `rc_calibration.json` | Numerical checks, comparison errors and fitted RC parameters |
| `interactive_report.html` | Standalone Plotly viewer, generated when Plotly is installed |
| `OpenChipCoolLoop.mo` | Exported Modelica loop |

The mesh-study command additionally produces `mesh_study_aligned.csv`,
`mesh_study_aligned.png` and `mesh_gci.json`. Large HTML and field files can
be generated locally with the commands above.

## OpenModelica

After generating `results/OpenChipCoolLoop.mo`, an installed `omc` can run
the model and compare its final die temperatures with the Python ODE:

```bash
python scripts/verify_openmodelica.py
```

Alternatively, `omc modelica/simulate.mos` runs the supplied simulation
script. The exported model is a lumped thermal network; hydraulic
`Modelica.Fluid` components are part of the roadmap. The verification
script returns exit code 2 when `omc` is unavailable.

## Microchannel cold plates

Two geometry profiles are available:

- **Literature benchmark:** 21 copper/water channels based on Qu & Mudawar
  (2002), with nominal channel dimensions of 231 × 713 µm.
- **36 mm concept:** 108 channels in a separate 36 × 36 mm plate at the
  reference conditions of 800 W, 30 °C inlet and 4 L/min.

![21-channel copper/water cross-section](results/cfd/geometry_paper.png)

![36 mm cold-plate cross-section](results/cfd/geometry_design36.png)

Generate a case and literature-reference plots using Python:

```bash
python scripts/prepare_openfoam.py --profile paper --level coarse --Re 890 --qeff 100 --case cfd_cases/qm890
python scripts/build_validation_matrix.py --output cfd_cases/validation_matrix
python scripts/cfd_reference_report.py
python scripts/plot_figure_references.py
```

The matrix contains six operating-point families at three mesh resolutions.
On a host with **OpenFOAM.com v2312** installed and its environment loaded,
the runner performs meshing, region checks and a steady CHT solve:

```bash
python scripts/run_openfoam_matrix.py cfd_cases/validation_matrix/qm_q100_Re890_coarse --execute
```

**Current CFD status:** generated case inputs and analytical/literature
references are available. OpenFOAM mesh checks, converged CHT fields and
experimental validation are pending; solver compatibility remains unverified.
The plots in `results/cfd/` contain geometry or analytical/literature data.

See the [CFD methodology](docs/qu_mudawar_openfoam_protocol.md),
[validation plan in Russian](docs/validation_plan_RU.md) and
[component status](docs/cfd_status.md) for details.

## Scope and limitations

- RC states represent average die temperatures. Local hotspot maxima
  come from the FEM field and can exceed the averages substantially.
- FEM training uses a fixed coolant boundary temperature. Dynamic coolant
  and low-flow scenarios extend beyond those training conditions;
  pump-fault predictions are qualitative.
- Package thermal properties are constant. Prepared CFD cases include
  tabulated water transport properties; OpenFOAM execution is pending.
- The CFD geometry covers the channel core. Manifolds, thermocouple holes,
  heater slots and cover conduction are not yet represented in full.
- Numerical verification and FEM/RC agreement do not establish experimental
  validation. The 95 °C threshold is an illustrative scenario parameter.

See [sources and verification](docs/literature_and_vv.md) for model provenance.

## Roadmap

- Run and verify the OpenFOAM benchmark across the prepared mesh levels.
- Extend the benchmark geometry with manifolds, cover conduction and the
  experimental thermocouple and heater details.
- Compare pressure and temperature predictions with independently digitized
  published observations and quantify uncertainty.
- Derive heat-transfer and pressure-drop curves for the 36 mm plate and
  couple them to the package and cooling-loop models.
- Cross-check the thermal FEM against an independent solver.
- Extend the Modelica loop with pump, pipe, valve and heat-exchanger hydraulics.
- Add contact-property uncertainty and chip-to-package thermomechanics.

## Documentation

- [Обзор проекта на русском](docs/summary_RU.md)
- [Sources, model provenance and numerical verification](docs/literature_and_vv.md)
- [OpenFOAM benchmark methodology](docs/qu_mudawar_openfoam_protocol.md)
- [План экспериментальной валидации](docs/validation_plan_RU.md)
- [CFD component status](docs/cfd_status.md)

## Scientific references

- Qu & Mudawar (2002), microchannel pressure drop and heat transfer:
  [DOI](https://doi.org/10.1016/S0017-9310(01)00337-4).
- Lee, Garimella & Liu (2005), water microchannel experiments:
  [DOI](https://doi.org/10.1016/j.ijheatmasstransfer.2004.11.019).
- Lee & Garimella (2006), thermally developing flow correlations:
  [DOI](https://doi.org/10.1016/j.ijheatmasstransfer.2006.02.011).
- Celik et al. (2008), discretization uncertainty and GCI:
  [DOI](https://doi.org/10.1115/1.2960953).

## License

[MIT](LICENSE). Dependencies retain their respective licenses.
