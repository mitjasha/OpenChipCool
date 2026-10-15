# Qu & Mudawar benchmark: validation methodology

The OpenFOAM benchmark targets the microchannel experiment of Qu & Mudawar
(2002). Case inputs and reference-data processing are implemented; meshing,
CHT execution and experimental validation are pending.

## Sources and geometry

Qu, W. & Mudawar, I. (2002), *Experimental and numerical study of pressure drop
and heat transfer in a single-phase micro-channel heat sink*, IJHMT **45**,
2549–2565. [DOI](https://doi.org/10.1016/S0017-9310(01)00337-4),
[author PDF](https://engineering.purdue.edu/mudawar/files/articles-all/2002/2002_1.pdf).

The experiment and the paper's numerical unit cell have different geometry.
The generated OpenChipCool case represents the full-length parallel-channel
core. Experimental manifolds, pressure taps, heater slots, thermocouple holes
and the conducting cover are outside its current scope.

| Parameter | Published reference | Generated core |
|---|---|---|
| Channels | 21, Fig. 2 | 21 water passages |
| Channel width/depth | 231 × 713 µm | Same nominal dimensions |
| Channel pitch | 467 µm, Table 1 | Same; remaining width assigned to outer edges |
| Copper footprint | 10 × 44.8 mm, Fig. 2; length 44.764 mm, Table 1 | 10 × 44.764 mm |
| Copper below channels | 5637 µm, Table 1 | Same |
| Inlet temperature | 15 °C | 288.15 K |
| Heat flux | 100/200 W/cm² of planform area | Equal total power applied at the copper base |
| Pressure observations | Between deep inlet/outlet plenums | Channel-core pressure only |
| Cover and heater slots | Included in the experimental fixture | Omitted |
| Water transport | Temperature-dependent viscosity in published CFD | IAPWS-based tabulated viscosity/conductivity prepared in v0.4; runtime compatibility pending |
| Copper thermocouples | Four positions, Table 1 and Fig. 7 | Approximate probe coordinates; transverse position and junction depth remain unverified |

The stated experimental uncertainties include wattmeter <0.5%, flowmeter <4%,
pressure <3.5% and thermocouple <0.3 °C. CSV references in `data/literature/`
are approximate readings of printed curves, with separate digitization
uncertainty. They are not the authors' raw measurements. The pressure CSV uses
an indicative ±0.05 bar reading uncertainty. Fig. 7 visual estimates are
available; the separate `tc_digitization_PENDING.csv` records the pending
higher-resolution digitization.

## Case generation

From an installed Python checkout:

```bash
python scripts/prepare_openfoam.py --profile paper --level coarse --Re 890 --qeff 100 --case cfd_cases/qm_890_coarse
python scripts/build_validation_matrix.py --output cfd_cases/validation_matrix
```

Outputs include `blockMeshDict`, region/property dictionaries, `geometry.json`,
`mesh_audit.json` and `setup_status.json`. The 21-channel core contains 86 blocks:
43 lower copper blocks and 43 channel-height blocks. Copper and water share
63 block faces. Python topology checks cover outer-face labels, shared
interfaces and non-manifold connectivity; OpenFOAM mesh quality is a separate
execution step.

| Profile | Coarse cells | Medium cells | Fine cells |
|---|---:|---:|---:|
| 21-channel benchmark | 21,400 | 171,200 | 1,369,600 |
| 108-channel, 36 mm plate | 108,400 | 867,200 | 6,937,600 |

These are specified cell counts. Each refinement doubles the subdivisions in
all three directions. The validation matrix has six operating-point families
and three mesh levels per family.

## OpenFOAM execution

The configuration targets **OpenFOAM.com v2312** and the steady solver
`chtMultiRegionSimpleFoam`. With the OpenFOAM environment loaded:

```bash
python scripts/run_openfoam_matrix.py cfd_cases/validation_matrix/qm_q100_Re890_coarse --execute
```

The runner calls `blockMesh`, `splitMeshRegions -cellZonesOnly -overwrite`,
regional `checkMesh` and the solver, retaining logs in the case directory.
Zone-based splitting groups the separate water ducts into one fluid region.
Missing tools or failed checks produce a nonzero exit code. A completed solver
process is recorded separately from convergence and validation.

Upstream documentation: [heat-transfer solvers](https://doc.openfoam.com/2312/tools/processing/solvers/rtm/heat-transfer/),
[tabulated transport](https://doc.openfoam.com/2312/tools/processing/models/thermophysical/transport/rtm/tabulated/).
Solver dictionaries and coupled boundary conditions have not yet been verified
by an OpenFOAM run.

## Observables and comparison

- **Fig. 4, pressure:** core pressure is the difference of area-averaged inlet
  and outlet pressures. Experimental tap-to-tap pressure additionally contains
  inlet/outlet plenum losses. `pressure_budget` implements the component
  relations in Eqs. (19)–(22); independent areas and velocities define those
  corrections. Calibration on the comparison points would remove their
  independence as validation data.
- **Fig. 5, coolant heating:** the observable is mass-flow-weighted outlet
  temperature minus inlet temperature. The associated energy check is
  `Q_water = mdot cp (Tout - Tin)`.
- **Fig. 7, copper temperatures:** the observables are four copper thermocouple
  rises above inlet temperature. Approximate x positions are 5, 16.588,
  28.176 and 39.764 mm. The prepared z coordinate of 3.175 mm and transverse
  locations remain subject to geometric verification. Channel-wall temperature
  and embedded thermocouple temperature are different observables.

The numerical study covers residual histories, mass and energy closure and
three systematically refined meshes. GCI is meaningful for monotone data with
a positive observed order and an established asymptotic regime. Three meshes
alone do not establish that regime.

Current study targets are mass imbalance <0.1%, energy imbalance <1% and
GCI <5%; comparison targets are mean pressure discrepancy <10% and average
absolute thermocouple discrepancy <5 K, with measurement and digitization
uncertainty included in the interpretation. These are project study targets,
not published Qu & Mudawar acceptance limits, and have not been achieved yet.

`scripts/extract_openfoam_result.py` reads solver postprocessing, checks mass
and energy closure and exports CSV. Its required arguments include the case,
`--plenum-correction-pa`, `--plenum-justification` and `--out`. The correction
and its justification come from an independent calculation for the fixture.
`scripts/validate_all_figures.py --cfd-export <combined.csv>` then produces
the comparison report. The legacy pressure-only interface is
`scripts/compare_cfd_experiment.py --cfd-results <export.csv>`.

## Analytical references

```bash
python scripts/cfd_reference_report.py
python scripts/plot_figure_references.py
```

The pressure plot uses a constant-viscosity, fully developed rectangular-duct
relation. Temperature plots contain literature readings and an energy-balance
reference. These outputs contain no OpenFOAM solution. Differences between the
channel-only pressure estimate and measured total pressure remain unresolved;
temperature-dependent viscosity and omitted manifolds are relevant model gaps.

## Transfer to the 36 mm plate

```bash
python scripts/prepare_openfoam.py --profile design36 --level coarse --case cfd_cases/plate36
```

The separate plate has 108 channels of 231 × 713 µm, 332 µm pitch, a 3 mm
copper base and a 36 × 36 mm footprint. Reference conditions are 800 W,
30 °C inlet and 4 L/min, corresponding to inlet Re ≈1632. Its geometry and
operating conditions differ from the published experiment.

The planned coupling uses geometry-specific pressure-drop curves, heat to
water, cold-plate thermal resistance and an equivalent planform heat-transfer
coefficient. Their definitions include reference coolant temperature, spatial
averaging and uncertainty. These curves will replace the current literature
surrogate after independent numerical checks and comparison.

See [CFD status](cfd_status.md), the [Russian validation plan](validation_plan_RU.md)
and [model provenance](literature_and_vv.md).
