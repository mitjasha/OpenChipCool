# Numerical verification

OpenChipCool separates numerical checks from experimental comparison.

- FEM checks cover source normalization, stationary and transient energy
  balances and geometry-aligned mesh refinement.
- The RC network is fitted to training transients and evaluated on a separate
  power sequence. Its states represent die averages; FEM peaks are reported
  separately.
- GCI estimates are conditional on monotone refinement and a positive observed
  order; three meshes alone do not confirm an asymptotic regime.
- CFD Python checks cover generated geometry and case inputs. OpenFOAM mesh
  checks, converged CHT solutions and experimental validation remain pending.
- OpenModelica execution checks require an installed `omc` and compare compiled
  outputs against the Python reference.

The [sources and model-provenance document](literature_and_vv.md) describes
these checks and their interpretation. The [CFD methodology](qu_mudawar_openfoam_protocol.md)
defines the pressure and temperature observables for the planned literature
comparison.
