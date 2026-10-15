# OpenChipCool v0.2: Sources, model provenance and V&V plan

## Source 1: direct experimental measurements + CFD (primary)

Qu, W.; Mudawar, I. (2002). *Experimental and numerical study of pressure drop
and heat transfer in a single-phase micro-channel heat sink*. International
Journal of Heat and Mass Transfer **45**, 2549–2565.
DOI: https://doi.org/10.1016/S0017-9310(01)00337-4
Open author PDF: https://engineering.purdue.edu/mudawar/files/articles-all/2002/2002_1.pdf

**Experimental**: oxygen-free copper coldplate; 21 channels, width 231 um,
height 713 um, 44.8 mm long, 10 mm wide footprint; water at approx 15 °C,
100 and 200 W/cm2 planform heat flux; measured pressure loss (Fig. 4a),
fluid outlet temperature (Fig. 5), thermocouple temperatures (Fig. 6).
139 <= Re <= 1672 in the 100 W/cm2 experiments; 385 <= Re <= 1289
for 200 W/cm2. Fig 4 experiment readings are **visually approximated** in
`data/literature/qu_mudawar_2002_pressure_approx.csv`, NOT original raw data
and NOT a hydraulics curve for the new 36-mm plate. Errors about 0.05 bar
are only digitization estimates, not source-supplied instrument uncertainties.

**CFD**: Figure 11(c) is the AUTHORS' conjugate-CFD output, not measurement.
It shows axially averaged Nu on bottom and side walls at Re=890, 100 W/cm2.
Our visually approximated traces in `qu_mudawar_2002_nu_cfd_fig11_approx.csv`
carry rough *reading* uncertainty of 20% and require independent digitization
for serious usage. The coldplate surrogate integrates these axial Nu profiles.
At other Re, Graetz-like `Re^(1/3)` scaling of entry enhancement is **our own
assumption**, not a Qu/Mudawar fitted correlation. At Re below 139 we apply a
linear reduction to zero. This is **unvalidated model-form extrapolation**.

The proposed 36 x 36 mm plate with 108 channels is **not the published test
article**: it shares nominal w/d geometry but differs in n, length, coolant
inlet and total power. Our equivalent planform coefficient assumes ideal fins
and a uniform fixed bulk-fluid temperature:

    h_equiv = (k_water / D_h) * [Nu_bottom*A_channel_bottom
                                    + Nu_side*A_channel_sides] / A_planform

It should NOT be called a measured heat-transfer coefficient for our package.
It is a literature-CFD-anchored *idealized Robin equivalent*. Wall-dependent
film coefficients, copper fin efficiency, coolant heating and manifolds must
be modeled before making physical design claims. Experimental pressure losses
from Fig 4 are stored *for independent benchmarking only*, not transferred
as the pressure drop of a 108-channel system.

## Source 2: experimental cross-check

Lee, P.-S.; Garimella, S. V.; Liu, D. (2005). *Investigation of heat
transfer in rectangular microchannels*. International Journal of Heat and
Mass Transfer **48**, 1688–1704.
DOI: https://doi.org/10.1016/j.ijheatmasstransfer.2004.11.019
Water experiments, 194–534 um widths, Re about 300–3500. Numerics vs measured
temperature trends average about 5% deviation, with inlet-condition matching.
Do not transfer their results as exact validation of our new geometry.

## Source 3: thermally developing CFD correlation research

Lee, P.-S.; Garimella, S. V. (2006). *Thermally developing flow and heat
transfer in rectangular microchannels of different aspect ratios*.
International Journal of Heat and Mass Transfer **49**, 3060–3067.
DOI: https://doi.org/10.1016/j.ijheatmasstransfer.2006.02.011
A better future model should implement paper's full Nu correlation with the
same boundary condition (uniform wall T vs uniform q'') and channel aspect
ratio, NOT interpret Re^(1/3) as a fitted expression from this source.

## Source 4: discretization verification standards

Celik, I. B. et al. (2008), *Procedure for Estimation and Reporting of
Uncertainty Due to Discretization in CFD Applications*, ASME J Fluids Eng
**130** 078001. DOI: https://doi.org/10.1115/1.2960953

Roache, P. J. (1994), *Perspective: A Method for Uniform Reporting of Grid
Refinement Studies*, ASME J Fluids Eng **116**, 405–413.
DOI: https://doi.org/10.1115/1.2910291

`mesh-study --refinements 1,2,3` refines all spatial axes and aligns
copper/TIM/silicon interfaces and exact chip edges. It computes (1) a
surface-quadrature die temperature, (2) maximum mesh-node top temperature,
(3) energy conservation and correct silicon source volume. Classical
Richardson + safety-factor 1.25 GCI are calculated only for monotone cases
with positive observed order. For noisy or oscillatory data GCI is flagged
inconclusive. It is NOT a guarantee of asymptotic behavior with three grids,
and peak values require extra care. Repeat with denser meshes / independent
solver before reporting a certified temperature.

## OpenModelica status and execution

The generator now exports the same CFD-Nu-anchored h(flow) dependence as the
Python lumped model, including the *explicitly unvalidated* low-flow regime.

```
python -m openchipcool.cli benchmark --refinement 2
python scripts/verify_openmodelica.py
```

`verify_openmodelica.py` checks for `omc`, calls the official scripting
API `loadFile/checkModel/simulate`, checks errors and requires a time-series
CSV output. If `omc` isn't installed it exits code 2 (UNVERIFIED, not PASS).
Install from https://openmodelica.org/download/download-linux/ and re-run.

**Modelica runs a 9-state lumped ODE**, not Modelica.Fluid pipes; it contains
no hydraulic pump pressure curve or inlet/outlet water transport. Numerical
comparison against Python should be done after confirming compiled Modelica
simulation and aligning time integrators / power scenarios. Static text tests
of a `.mo` file never certify compilation or simulation.
