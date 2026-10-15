# Qu & Mudawar (2002) — 21-channel OpenFOAM CHT reproduction protocol

**OpenChipCool v0.3 status: GEOMETRY GENERATED; CFD NOT RUN; NOT EXPERIMENTALLY VALIDATED.**

## Publication and source traceability

Qu, W. & Mudawar, I. (2002), *Experimental and numerical study of pressure drop and heat transfer in a single-phase micro-channel heat sink*, International Journal of Heat and Mass Transfer **45**, 2549–2565.
DOI: https://doi.org/10.1016/S0017-9310(01)00337-4
Full article: https://engineering.purdue.edu/mudawar/files/articles-all/2002/2002_1.pdf

The paper's experimental fixture (Fig. 2) is not identical to its unit-cell CFD model (Fig. 3/Table 1). Our first-stage mesh reproduces the **21 full-length parallel channel core**, NOT the experimental inlet/outlet plenums, heater-body transverse slots, pressure taps, or conducting polycarbonate cover. Its core-only pressure drop must NOT be compared directly against the published measured plenum-to-plenum pressure as if they were the same quantity.

| Parameter | Published, citation | Present model |
|---|---|---|
| Channels | 21, Fig. 2 | 21 fluid passages |
| Channel size | 231 µm wide × 713 µm deep, abstract/Fig. 2 | same |
| Unit-cell channel pitch | 118 + 231 + 118 = 467 µm, Table 1 | same; remainder of 10-mm plate distributed at two edges |
| Copper footprint | 10 × 44.8 mm (Fig. 2); L=44.764 mm (Table 1) | 10 × 44.764 mm |
| Copper below channels | 5637 µm (Table 1) | same |
| Water inlet | 15°C, paper p.2553 | 288.15 K |
| Applied heat flux | 100 and 200 W/cm² relative to TOP PLANFORM, Eq. 1 | fixed heat flux on copper base; 447.64 and 895.28 W |
| Experimental flow | 93.5–1136.4 mL/min (100), 259.1–859.5 mL/min (200), p.2553 | corresponding Re via Eq. 2–3 |
| Measured Reynolds | 139–1672 (100), 385–1289 (200), p.2553 | operating points inside these ranges |
| Pressure measurement | between inlet and outlet DEEP PLENUMS, p.2553, Eq. 18–22 | **only core channels**; missing manifold losses |
| Polycarbonate cover / heater slots | Fig. 2, upper convection h=10 W/m²K in CFD, p.2555 | adiabatic top, simplified bottom heat flux; not exact experiment |
| Water viscosity | numerical unit cell uses temperature-varying water viscosity, p.2554 | constant liquid viscosity at inlet; can bias pressure predictions |
| Four thermocouples | Fig. 2 and Table 1 | suggested x = 5, 16.588, 28.176, 39.764 mm; depth requires checking original thermocouple plane |

Paper's stated experimental errors: wattmeter <0.5%, flowmeter <4%, pressure <3.5%, thermocouple <0.3°C. Largest explicitly stated **measured** pressure: 0.86 bar, Re=1672, q''=100 W/cm², p.2556. **Our digitized pressure-marker CSV is an approximate hand reading of Fig. 4(a), plus more uncertain Fig. 4(b) markers; it is NOT the authors' original raw data.** Maintain a separate ±0.05-bar *digitization* error assumption. Temp measurements in Fig. 7 have NOT yet been extracted; do not substitute guessed temperature numbers.

Additional scientific reference: Celik et al. (2008), *Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications*, Journal of Fluids Engineering 130, 078001. https://doi.org/10.1115/1.2960953 (GCI methodology, not proof of numerical convergence).

OpenFOAM upstream docs: https://doc.openfoam.com/2312/tools/processing/solvers/rtm/heat-transfer/chtMultiRegionFoam/ and https://doc.openfoam.com/2606/tools/pre-processing/mesh/generation/blockMesh/ . Config targets **OpenFOAM.com v2312 steady-state `chtMultiRegionSimpleFoam`**; its solver compatibility has not been tested here.

## Reproducible geometry / case generation

Install Python requirements from repository root:

```bash
pip install -e '.[dev]'
pytest -q
python scripts/prepare_openfoam.py --profile paper --level coarse --Re 890 --qeff 100 --case cfd_cases/qm_890_coarse
```

This creates `system/blockMeshDict`, `constant/regionProperties`, per-region material and solver dictionaries, `geometry.json`, `mesh_audit.json` and `setup_status.json`. The block mesh is **86 blocks**: 43 bottom copper blocks and 43 channel-height blocks, consisting of copper fins interleaved with 21 water passages. Copper and water share hexahedral block faces. The combinatorial audit checks that all outer faces are labeled exactly once, no non-manifold block interface exists, and 63 water/copper block-faces are shared.

On a host where OpenFOAM.com tools are actually installed and its environment sourced:

```bash
python scripts/prepare_openfoam.py --profile paper --level coarse --Re 890 --qeff 100 --case cfd_cases/qm_890_coarse --mesh
cd cfd_cases/qm_890_coarse
chtMultiRegionSimpleFoam > log.cht 2>&1
```

`--mesh` attempts **blockMesh**, **splitMeshRegions -cellZonesOnly -overwrite**, `checkMesh` for water and copper and initialization of actual split-region boundary fields. We intentionally use `-cellZonesOnly`: ordinary connected-component splitting separates the 21 unconnected water ducts into 21 separate regions. Solver dictionaries and coupled interface boundary condition names must still be tested with the installed OpenFOAM build. A solver exit status alone does NOT establish steady convergence or experimental validity. If `blockMesh`, `splitMeshRegions`, `checkMesh` or solver is missing, the script returns an explicit error rather than reporting PASS.

### Grid levels

| Case | Coarse | Medium | Fine |
|---|---:|---:|---:|
| 21-channel benchmark, combined hexes | 21,400 | 171,200 | 1,369,600 |
| 108-channel 36-mm concept, combined hexes | 108,400 | 867,200 | 6,937,600 |

These are **planned cell counts**, not logs from `checkMesh`. Each refinement doubles x/y/z subdivisions. Solve at least 3 levels, check actual face matching/non-orthogonality and report GCI of pressure and selected thermocouple temperatures where convergence is monotone and asymptotic. For example use operating cases: q''=100 W/cm² at Re 226/442/890/1454/1672, and q''=200 W/cm² at Re 442/685/864/1220. Use exact flow rates computed in `experimental_condition.json`; align to Fig. 4/5/7 test conditions when digitizing.

## Physics, measurements, and validity criteria

- Fluid: incompressible laminar flow and convective-diffusive energy equation; copper: steady conduction; coupled interface conserves `T` and normal heat flux. No phase change.
- Inlet mass flow is set using Qu & Mudawar Eq. (2–3), `Re = rho u Dh/mu`, `u = Vdot / (N A_channel)`. Copper base carries specified heat flux on planform area.
- Track `Delta p_core = <p>_channel,in - <p>_channel,out` separately from the experiment's **tap-to-tap** `Delta p_total`. Model or independently determine entrance/exit manifold losses per Eq. (18–22), or add the physical manifolds to the mesh; do **not** equate the two.
- Obtain coolant `Tout` with **mass-flow-weighted** outlet averaging, not unweighted face mean. Check `Qfluid = mdot cp (Tout - Tin)` against `Qheater`; compare experimental outlet rise from Fig. 5. Sample four copper **thermocouple** coordinates for comparison with Figs. 6–7 (Fig. 7 shows `T_thermocouple − T_inlet`), not the channel surface.
- Target numerical verification criteria (proposed, NOT passed): total mass imbalance <0.1%; energy imbalance <1%; GCI <5% for pressure/temperature in asymptotic regime. Proposed correlation objectives: plenum-corrected mean pressure discrepancy <10% and average absolute thermocouple error <5°C, after accounting for stated uncertainties. These targets do not originate in Qu & Mudawar and are not safety acceptance limits.

For a **theoretical** check without OpenFOAM:

```bash
python scripts/cfd_reference_report.py
```

The generated `results/cfd/pressure_reference_analytical_only.png` compares hand-digitized **paper** pressure markers to a **core-only constant-viscosity analytical rectangular-duct relation**. The black curve is NOT OpenFOAM, and a discrepancy (even a channel-only prediction above the measured full Δp) must be reported as such. Temperature-dependent viscosity and the omitted plenum require investigation; do not fit fictitious local losses to disguise discrepancies.

Export actual *independent solver outputs* as columns `Re,qeff_W_cm2,dp_core_Pa,dp_total_Pa,Tout_C,energy_error_fraction`, with the corrected Δp's provenance documented. Then compare:

```bash
python scripts/compare_cfd_experiment.py --cfd-results path/to/your_actual_OpenFOAM_export.csv
```

This command fails if the CFD export is absent. Fig. 7 copper thermocouple observations still need explicit documented graphical digitization (e.g. WebPlotDigitizer) before reporting numerical temperature errors. A successful pressure plot by itself is NOT validation of conjugate heat transfer.

## Reusing the methodology on the new OpenChipCool cold plate

```bash
python scripts/prepare_openfoam.py --profile design36 --level coarse --case cfd_cases/plate36 --mesh
```

**Design assumptions:** 36 × 36-mm planform, 108 channels of 231 × 713 µm, 332 µm pitch, 3-mm copper below slots, water inlet 30°C, power 800 W = 61.7284 W/cm², flow 4 L/min (inlet Re ≈ 1632 at prescribed 30°C properties). These are **not** Qu & Mudawar dimensions or operating points and cannot be labeled experimentally validated based on the original 21-channel study. This is a *channel-core* model without manifolds or real four-chiplet spatially varying heat flux. It does not yet define the full 36-mm cold-plate hydraulic design.

After obtaining converged and reviewed CHT solutions **for the new geometry**, export `Delta p(flow)` (accounting for manifolds), energy-conserving `Q_to_water`, `Rth_coldplate(flow) = (Tbase_avg-Tcoolant_ref)/Qfluid` and area-averaged `h_equiv = Qfluid/[Aplanform*(Tbase_avg-Tcoolant_bulk)]`. Report the reference temperature and spatial averaging procedure explicitly. These characteristics can then replace the placeholder `h(flow)` in the FEM/RC/Modelica loop, ideally accompanied by a fitted uncertainty band and heat-flux maps. **Not connected yet: no CFD solution exists in the archive.**

Overall intended evidence chain: **published experiment → generated 21-channel case → actual OpenFOAM mesh/CHT + GCI → comparison with Δp and four copper temperatures → 36-mm CHT → FEM package heat flux → RC model → OpenModelica cooling loop.** The chain currently stops after geometry generation and analytical screening.
