"""Literature-anchored coldplate surrogate (NOT own CFD or experimental data).

Qu & Mudawar, IJHMT 45 (2002) 2549-2565, Fig. 11(c), DOI:
https://doi.org/10.1016/S0017-9310(01)00337-4
Digitized *approximate* CFD Nu at Re=890 (not measurements), channel bottom
and side walls. Other Reynolds numbers use an assumed Graetz scaling for the
entrance enhancement and are NOT directly confirmed by the article. Flow
below Re=139 uses an explicit shutdown extrapolation (not validated).

Geometric scaling to a *hypothetical* 36-mm plate with 108 parallel channels
is a new design assumption, NOT a result from Qu & Mudawar. A flat Robin h
is only an idealized equivalent coefficient neglecting fin-efficiency,
axial bulk warming, plenum maldistribution and TIM/contact resistance.
"""
from dataclasses import dataclass
import numpy as np

PAPER = 'Qu & Mudawar (2002), doi:10.1016/S0017-9310(01)00337-4, Fig. 11(c)'
# Manually read from Fig. 11(c), approximate, ~20% reading uncertainty; not
# directly supplied as a table in the publication. x = 0 ... 44.8 mm.
X = np.array([0, .002, .005, .010, .020, .030, .0448])
NU_BOTTOM = np.array([19.0, 12.0, 8.1, 6.3, 5.0, 4.6, 4.3])
NU_SIDE = np.array([30.0, 17.0, 11.0, 8.7, 7.8, 7.7, 7.4])

@dataclass(frozen=True)
class ChannelDesign:
    n_channels: int = 108  # hypothetical 36 mm plate; NOT Qu/Mudawar device
    width_m: float = 231e-6  # adopted from experimental geometry
    depth_m: float = 713e-6  # adopted from experimental geometry
    length_m: float = .036   # new model geometry
    plate_width_m: float = .036
    plate_length_m: float = .036
    water_density: float = 995.7  # kg/m^3 approx 30 degC: engineering property assumption
    water_viscosity: float = 0.000798  # Pa*s approx 30 degC
    water_conductivity: float = .613  # W/m/K approx 30 degC
    reference_re: float = 890.  # study's Fig. 11 CFD condition
    min_measured_re: float = 139.  # Qu & Mudawar experimental range

    @property
    def hydraulic_diameter(self):
        return 2*self.width_m*self.depth_m/(self.width_m+self.depth_m)

    @property
    def planform_area(self):
        return self.plate_width_m*self.plate_length_m

    @property
    def bottom_area(self):
        return self.n_channels*self.width_m*self.length_m

    @property
    def side_area(self):
        return 2*self.n_channels*self.depth_m*self.length_m


def reynolds(flow_lpm: float, plate=ChannelDesign()):
    if flow_lpm < 0:
        raise ValueError('Negative coldplate flow')
    m3s = flow_lpm*1e-3/60
    velocity = m3s/(plate.n_channels*plate.width_m*plate.depth_m)
    return plate.water_density*velocity*plate.hydraulic_diameter/plate.water_viscosity


def reference_nusselt():
    # Integrate axial profile rather than use an arbitrary power-law h.
    return (float(np.trapezoid(NU_BOTTOM,X)/(X[-1]-X[0])),
            float(np.trapezoid(NU_SIDE,X)/(X[-1]-X[0])))


def equivalent_htc(flow_lpm: float, plate=ChannelDesign()):
    """Equivalent planform-area HTC from approximate CFD Nu and geometry.

    Under 139 Re, suppress heat transfer to zero with *explicitly unvalidated*
    linear shutdown; for Re above CFD anchor extrapolate entrance Nu ~ Re^1/3.
    Assumes perfect fin conduction and uniform bulk coolant temperature.
    """
    re = reynolds(flow_lpm,plate)
    if re == 0:
        return 0.
    b,s = reference_nusselt()
    # Asymptotic values estimated from downstream CFD curves, NOT a published
    # formula: variation beyond the CFD reference condition is an assumption.
    graetz = (re/plate.reference_re)**(1./3.)
    nu_bottom = 3.9+(b-3.9)*graetz
    nu_side = 6.8+(s-6.8)*graetz
    k_over_d = plate.water_conductivity/plate.hydraulic_diameter
    h = k_over_d*(nu_bottom*plate.bottom_area+nu_side*plate.side_area)/plate.planform_area
    # In the pump-fault regime the constant-bulk forced convection hypothesis
    # is invalid. Explicit model-form uncertainty, not fitted to data.
    if re < plate.min_measured_re:
        h *= re/plate.min_measured_re
    return float(h)


def heat_transfer_coefficient(flow_lpm: float, spec):
    # Physical package area matches our hypothetical 36-mm coldplate planform.
    plate=ChannelDesign(plate_width_m=spec.side,plate_length_m=spec.side,
                        length_m=spec.side)
    return equivalent_htc(flow_lpm,plate)


def summary(flow_lpm=4.,plate=ChannelDesign()):
    b,s=reference_nusselt()
    return dict(source=PAPER, digitalization='visual approximate, not original raw data',
                flow_lpm=flow_lpm, reynolds=reynolds(flow_lpm,plate),
                nu_bottom_ref=b,nu_side_ref=s, h_equiv_W_m2K=equivalent_htc(flow_lpm,plate),
                validity=('within experimental Re range' if 139<=reynolds(flow_lpm,plate)<=1672
                          else 'EXTRAPOLATED: outside Qu/Mudawar Re range'),
                notes=['Nu curves are CFD predictions, not measured h',
                       '108 channels and 36-mm plate are proposed, not tested',
                       'Flow-dependence uses assumed Graetz entrance scaling',
                       'Unknown shutdown regime below Re=139',
                       'Perfect fin efficiency and fixed uniform coolant bulk; not package validation'])
