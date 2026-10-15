"""Experimental unexecuted OpenFOAM.com v2312 CHT setup for generated blockMesh.

Run after `blockMesh` and `splitMeshRegions -cellZonesOnly -overwrite`.
Keep explicit UNVERIFIED marker until a solver log, mesh check and validation
comparison are available. Other OpenFOAM distributions need adapted dictionaries.
"""
from __future__ import annotations
import re
from pathlib import Path
from .reference import experimental_condition
from .geometry import Geometry, geometry_profile


def _header(obj,klass='dictionary'):
    return f'FoamFile {{ version 2.0; format ascii; class {klass}; object {obj}; }}\n'


def _field(name, dimension, inner, patches, vector=False):
    a=_header(name,'volVectorField' if vector else 'volScalarField')
    return a+f'dimensions {dimension};\ninternalField uniform {inner};\nboundaryField\n{{\n'+''.join(
        f' {n}\n {{\n  {val}\n }}\n' for n,val in patches.items())+'}\n'


def _names(boundary:Path)->list[str]:
    if not boundary.is_file():
        raise FileNotFoundError(f'Missing split mesh: {boundary}; run blockMesh and splitMeshRegions first')
    s=boundary.read_text()
    names=re.findall(r'(?m)^\s*([A-Za-z_][\w]*)\s*\n\s*\{\s*\n\s*type\s+(?:\w+)',s)
    if not names:
        # OpenFOAM file format allows patch declaration/brace on a single line.
        names=re.findall(r'(?m)^\s*([A-Za-z_][\w]*)\s*\{\s*type\s+',s)
    return list(dict.fromkeys(n for n in names if n!='FoamFile'))


def _find_interface(name:str):
    return ('_to_' in name) or ('to_' in name) or name.startswith('interface')


def write_fv_dictionaries(case:Path, g:Geometry, re_num:float, flux:float, solver='chtMultiRegionSimpleFoam'):
    root=Path(case)
    condition=experimental_condition(re_num,flux,g)
    nx=g.nx
    ctl=_header('controlDict')+f'''application {solver};
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 2500;
deltaT 1;
writeControl timeStep;
writeInterval 200;
purgeWrite 2;
writeFormat ascii;
writePrecision 8;
runTimeModifiable yes;
'''
    # OpenFOAM.com v2312 style function objects. These are *requested*
    # diagnostics; they become evidence only when postProcessing files exist.
    tc_x=[.005,.016588,.028176,.039764] if g.name=='qu_mudawar_2002_core' else []
    probe_x='\n'.join(f'            ({x:.8g} {g.width_m/2:.9g} 0.003175)' for x in tc_x)
    # T at copper thermocouple x stations; heat sink slots are NOT resolved.
    funcs=f'''functions
{{
    tc_copper
    {{
        type probes;
        libs (sampling);
        region copper;
        fields (T);
        probeLocations
        (
{probe_x}
        );
        writeControl timeStep;
        writeInterval 20;
    }}
    inlet_p_areaAverage
    {{
        type surfaceFieldValue;
        libs (fieldFunctionObjects);
        region water;
        regionType patch;
        name inlet;
        operation areaAverage;
        fields (p);
        writeControl timeStep;
        writeInterval 20;
    }}
    outlet_p_areaAverage
    {{
        type surfaceFieldValue;
        libs (fieldFunctionObjects);
        region water;
        regionType patch;
        name outlet;
        operation areaAverage;
        fields (p);
        writeControl timeStep;
        writeInterval 20;
    }}
    outlet_mass_flow
    {{
        type surfaceFieldValue;
        libs (fieldFunctionObjects);
        region water;
        regionType patch;
        name outlet;
        operation sum;
        fields (phi);
        writeControl timeStep;
        writeInterval 20;
    }}
    inlet_mass_flow
    {{
        type surfaceFieldValue;
        libs (fieldFunctionObjects);
        region water;
        regionType patch;
        name inlet;
        operation sum;
        fields (phi);
        writeControl timeStep;
        writeInterval 20;
    }}
    outlet_T_mass_weighted
    {{
        type surfaceFieldValue;
        libs (fieldFunctionObjects);
        region water;
        regionType patch;
        name outlet;
        operation weightedAverage;
        weightField phi;
        fields (T);
        writeControl timeStep;
        writeInterval 20;
    }}
}}
'''
    if not tc_x:
        # No Qu/Mudawar thermocouple arrangement on the new design concept.
        funcs=funcs[:funcs.index('    tc_copper')] + funcs[funcs.index('    inlet_p_areaAverage'):]
    (root/'system/controlDict').write_text(ctl+funcs)
    (root/'system').mkdir(exist_ok=True)
    for r in ('water','copper'):
        sd=root/'system'/r;sd.mkdir(exist_ok=True,parents=True)
        if r=='water':
            schemes='''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
 default Gauss upwind;
 div(phi,U) Gauss linearUpwind grad(U);
 div(phi,h) Gauss upwind;
 div(phi,e) Gauss upwind;
 div(phi,K) Gauss upwind;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; p_rgh; }
'''
            solution='''solvers
{
 "(p|p_rgh)" { solver GAMG; tolerance 1e-8; relTol 0.01; smoother GaussSeidel; }
 "(U|h|e|T)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
 "(U|h|e|T)Final" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }
}
SIMPLE { nNonOrthogonalCorrectors 0; momentumPredictor yes; consistent yes; }
relaxationFactors
{
 fields { p 0.3; p_rgh 0.3; }
 equations { U 0.5; h 0.7; e 0.7; T 0.7; }
}
'''
        else:
            schemes='''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
'''
            solution='''solvers
{
 "(h|e|T)" { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0.01; }
 "(h|e|T)Final" { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0; }
}
SIMPLE { nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { h 0.8; e 0.8; T 0.8; } }
'''
        (sd/'fvSchemes').write_text(_header('fvSchemes')+schemes)
        (sd/'fvSolution').write_text(_header('fvSolution')+solution)
    (root/'system/fvSolution').write_text(_header('fvSolution')+'SIMPLE { nNonOrthogonalCorrectors 0; }\n')
    (root/'system/fvSchemes').write_text(_header('fvSchemes')+'ddtSchemes { default steadyState; }\n')
    return condition


def write_thermophysics(case:Path,g:Geometry|None=None, variable_water_transport:bool=True):
    p=Path(case)
    newplate=(g is not None and g.name=='openchipcool_36mm_concept')
    for r in ('water','copper'):
        d=p/'constant'/r;d.mkdir(parents=True,exist_ok=True)
        if r=='water':
            mu=0.000798 if newplate else 0.001138
            Pr=5.44 if newplate else 8.07
            cp=4178 if newplate else 4180
            rho=995.7 if newplate else 999.1
            therm='''thermoType
{
 type heRhoThermo;
 mixture pureMixture;
 transport const;
 thermo hConst;
 equationOfState rhoConst;
 specie specie;
 energy sensibleEnthalpy;
}
mixture
{
 specie { nMoles 1; molWeight 18.015; }
 thermodynamics { Cp 4180; Hf 0; }
 transport { mu 0.001138; Pr 8.07; }
 equationOfState { rho 999.1; }
}
'''
            therm=(therm.replace('Cp 4180',f'Cp {cp}').replace('mu 0.001138',f'mu {mu}')
                   .replace('Pr 8.07',f'Pr {Pr}').replace('rho 999.1',f'rho {rho}'))
            if variable_water_transport:
                from .water_properties import openfoam_tabulated_transport
                therm=therm.replace('transport const;', 'transport tabulated;')
                therm=therm.replace(f'transport {{ mu {mu}; Pr {Pr}; }}',openfoam_tabulated_transport())
                if 'transport tabulated;' not in therm or 'mu\n' not in therm:
                    raise AssertionError('Failed to insert IAPWS water transport table')
            (d/'turbulenceProperties').write_text(_header('turbulenceProperties')+'simulationType laminar;\n')
            (d/'g').write_text(_header('g','uniformDimensionedVectorField')+
                'dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n')
        else:
            therm='''thermoType
{
 type heSolidThermo;
 mixture pureMixture;
 transport constIso;
 thermo hConst;
 equationOfState rhoConst;
 specie specie;
 energy sensibleEnthalpy;
}
mixture
{
 specie { nMoles 1; molWeight 63.546; }
 thermodynamics { Cp 385; Hf 0; }
 transport { kappa 401; }
 equationOfState { rho 8960; }
}
'''
        (d/'thermophysicalProperties').write_text(_header('thermophysicalProperties')+therm)


def write_fields_after_split(case:Path, velocity:float, tin_K:float, flux_W_cm2:float):
    base=Path(case)
    water_names=_names(base/'constant/water/polyMesh/boundary')
    solid_names=_names(base/'constant/copper/polyMesh/boundary')
    if not {'inlet','outlet'}.issubset(water_names) or 'heater' not in solid_names:
        raise ValueError(f'Unexpected split region patches: water={water_names}, copper={solid_names}')
    if not any(_find_interface(n) for n in water_names) or not any(_find_interface(n) for n in solid_names):
        raise ValueError('No mapped solid/fluid interface in split meshes: check -cellZonesOnly option')
    p0=base/'0/water';p0.mkdir(parents=True,exist_ok=True)
    s0=base/'0/copper';s0.mkdir(parents=True,exist_ok=True)
    v=f'({velocity:.12g} 0 0)'
    def coupling(region):
        return 'type compressible::turbulentTemperatureCoupledBaffleMixed;\n  Tnbr T;\n  kappaMethod '+('fluidThermo' if region=='water' else 'solidThermo')+f';\n  value uniform {tin_K:.8g};'
    Up={}; Tp={}; pp={}; p_rgh={}; alphat={}
    for n in water_names:
        if n=='inlet':
            Up[n]=f'type fixedValue; value uniform {v};'
            Tp[n]=f'type fixedValue; value uniform {tin_K:.8g};'
            p_rgh[n]='type fixedFluxPressure; value uniform 0;'
        elif n=='outlet':
            Up[n]='type zeroGradient;'
            Tp[n]='type zeroGradient;'
            p_rgh[n]='type fixedValue; value uniform 0;'
        elif _find_interface(n):
            Up[n]='type noSlip;'
            Tp[n]=coupling('water')
            p_rgh[n]='type fixedFluxPressure; value uniform 0;'
        else:
            Up[n]='type noSlip;'
            Tp[n]='type zeroGradient;'
            p_rgh[n]='type fixedFluxPressure; value uniform 0;'
        pp[n]='type calculated; value uniform 0;'
        alphat[n]='type calculated; value uniform 0;'
    Ts={}
    for n in solid_names:
        if n=='heater':
            # +/- sign checked against heat flux outward: q_out=-k*dT/dn;
            # positive outward gradient -> imposed heat flows INTO copper.
            Ts[n]=f'type fixedGradient; gradient uniform {flux_W_cm2*1e4/401.:.12g};'
        elif _find_interface(n):
            Ts[n]=coupling('copper')
        else: Ts[n]='type zeroGradient;'
    fields={
        p0/'U':_field('U','[0 1 -1 0 0 0 0]','(0 0 0)',Up,vector=True),
        p0/'p':_field('p','[1 -1 -2 0 0 0 0]','0',pp),
        p0/'p_rgh':_field('p_rgh','[1 -1 -2 0 0 0 0]','0',p_rgh),
        p0/'T':_field('T','[0 0 0 1 0 0 0]',str(tin_K),Tp),
        p0/'alphat':_field('alphat','[1 -1 -1 0 0 0 0]','0',alphat),
        s0/'T':_field('T','[0 0 0 1 0 0 0]',str(tin_K),Ts),
    }
    for p,v in fields.items(): p.write_text(v)
    return sorted(map(str,fields))
