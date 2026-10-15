"""All geometry and material parameters are synthetic design assumptions (SI units)."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Spec:
    side: float = 0.036
    chip_side: float = 0.012
    chip_pitch: float = 0.014
    copper_thickness: float = 0.003
    tim_thickness: float = 0.0001
    silicon_thickness: float = 0.0005
    k_copper: float = 390.0
    k_tim: float = 5.0
    k_silicon: float = 130.0
    k_mold: float = 0.7
    volumetric_heat_copper: float = 3.45e6
    volumetric_heat_tim: float = 1.8e6
    volumetric_heat_silicon: float = 1.65e6
    volumetric_heat_mold: float = 1.4e6
    # h_ref is derived from Qu & Mudawar (2002) CFD Nu curves by the property below.
    flow_ref_lpm: float = 4.0
    # Below the reference experimental Re range, dry-out/pump-stop behavior is unknown.
    stagnant_fraction: float = 0.0
    coolant_initial: float = 30.0
    coolant_capacity: float = 42000.0  # J/K (~10 L water), effective lumped capacity
    radiator_UA: float = 300.0  # W/K, synthetic radiator/heat exchanger
    ambient: float = 25.0
    junction_limit: float = 95.0

    @property
    def h_ref(self):
        return h_from_flow(self.flow_ref_lpm, self)

    @property
    def z_tim(self):
        return self.copper_thickness

    @property
    def z_die(self):
        return self.copper_thickness + self.tim_thickness

    @property
    def height(self):
        return self.z_die + self.silicon_thickness


def chip_centers(spec: Spec):
    p = spec.chip_pitch / 2
    return [(-p, -p), (p, -p), (-p, p), (p, p)]


def h_from_flow(flow_lpm: float, spec: Spec):
    """CFD-Nusselt-anchored *effective* planform HTC, NOT coldplate validation.

    See literature.py for exact source and physically important limitations.
    """
    from .literature import heat_transfer_coefficient
    return heat_transfer_coefficient(flow_lpm, spec)
