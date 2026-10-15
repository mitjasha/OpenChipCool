"""OpenFOAM CHT case preparation and evidence-based CFD verification (v0.3)."""
from .geometry import Geometry, mesh_dictionary, geometry_metrics, geometry_profile
from .reference import laminar_core_pressure_pa, experimental_condition, compare_against_paper
