# CFD status

The CFD module currently provides case generation and reference-data processing.
OpenFOAM execution and experimental validation are pending.

| Component | Available material | Status |
|---|---|---|
| Literature benchmark | Qu & Mudawar (2002) geometry and measurement definitions | Documented |
| 21-channel core geometry | Parametric copper/water block geometry | Generated |
| Block topology | Python checks of shared faces, zones and external patches | Implemented |
| Validation matrix | Six operating-point families, three mesh levels each | Generated inputs |
| Water transport | IAPWS-based tabulated viscosity and thermal conductivity | Python implementation; OpenFOAM compatibility pending |
| OpenFOAM meshing | Meshing and regional checkMesh runner | Execution pending |
| Steady CHT solution | Solver dictionaries and run script | Execution pending |
| Pressure and temperature references | Approximate readings of published Figs. 4, 5 and 7 | Available with digitization uncertainty |
| CFD/experiment comparison | Export and comparison scripts | Awaiting solver outputs |
| CFD grid convergence | Three specified mesh resolutions | Awaiting solved meshes |
| 36 mm cold plate | 108-channel core geometry | Generated inputs |
| CFD coupling to FEM/RC/Modelica | Current loop uses a literature-based convection surrogate | Pending |

Geometry images and analytical/literature plots are available in `results/cfd/`.
They do not contain OpenFOAM solution fields. The channel-core geometry omits
parts of the experimental fixture, so a future core solve will have a narrower
scope than validation of the complete experiment.

See the [methodology](qu_mudawar_openfoam_protocol.md) and
[validation plan](validation_plan_RU.md) for observables and comparison criteria.
