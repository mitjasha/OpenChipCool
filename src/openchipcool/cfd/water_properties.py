"""IAPWS SR6-08(2011), liquid water at 0.1 MPa.

Eqs. (7) and (8) are published approximations of the IAPWS water viscosity
and conductivity models. Use only in stated ranges. Density and Cp still
constant in the CHT model and are independent modelling approximations.
DOI/source: https://iapws.org/technical-guidance/release/LiquidWater
"""
from __future__ import annotations
import math

IAPWS_URL='https://iapws.org/technical-guidance/release/LiquidWater'


def viscosity_pa_s(T_K:float)->float:
    if not (273.15<=T_K<=373.15):
        raise ValueError('This single-phase table excludes extrapolation outside 273.15–373.15 K')
    x=T_K/300
    return 1e-6*(280.68*x**(-1.9)+511.45*x**(-7.7)+61.131*x**(-19.6)+0.45903*x**(-40.0))


def conductivity_w_mk(T_K:float)->float:
    if not (273.15<=T_K<=373.15):
        raise ValueError('This single-phase table excludes extrapolation outside 273.15–373.15 K')
    x=T_K/300
    return 1.6630*x**(-1.15)-1.7781*x**(-3.4)+1.1567*x**(-6.0)-0.432115*x**(-7.6)


def openfoam_tabulated_transport()->str:
    temp=[273.15+5.*i for i in range(21)]
    mu='\n'.join(f'        ({t:.8f} {viscosity_pa_s(t):.12g})' for t in temp)
    k='\n'.join(f'        ({t:.8f} {conductivity_w_mk(t):.12g})' for t in temp)
    return f'''transport
{{
    mu
    (
{mu}
    );
    kappa
    (
{k}
    );
}}
'''
