"""Geometric/evidence checks: do not require or fake CFD solver availability."""
from pathlib import Path
import math
import pytest
from openchipcool.cfd.geometry import Geometry,geometry_profile,geometry_metrics,mesh_dictionary,write_geometry_case
from openchipcool.cfd.reference import experimental_condition,laminar_core_pressure_pa,compare_against_paper
from openchipcool.cfd.foam_setup import write_fv_dictionaries,write_thermophysics


def test_paper_geometry_and_mass_balance():
    g=geometry_profile('paper','coarse')
    m=geometry_metrics(g)
    assert m['channels']==21
    assert math.isclose(g.length_m,44.764e-3)
    assert math.isclose(g.channel_width_m,231e-6)
    assert m['water_volume_m3']>0
    assert m['cells_total']==m['cells_copper']+m['cells_fluid']
    # Exact geometric copper/fluid volume complement to bounding prism.
    assert math.isclose(m['water_volume_m3']+m['copper_volume_m3'],
                        g.width_m*g.length_m*(g.base_m+g.channel_depth_m),rel_tol=1e-13)


def test_openchipcool_design_36mm_and_distinct_from_experiment():
    g=geometry_profile('design36','coarse')
    m=geometry_metrics(g)
    assert g.width_m==g.length_m==.036
    assert g.n_channels==108
    assert m['cells_total']>geometry_metrics(geometry_profile('paper','coarse'))['cells_total']
    assert m['edge_margin_m']>0


def test_topology_foam_dictionary_all_faces_and_zones():
    g=geometry_profile('paper','coarse')
    s=mesh_dictionary(g)
    assert s.count('hex (')==2*(2*g.n_channels+1)
    assert s.count('simpleGrading')==2*(2*g.n_channels+1)
    assert 'inlet {' in s and 'outlet {' in s and 'heater {' in s
    assert ' copper (' in s and ' water (' in s
    assert s.count('(')>s.count('hex')


def test_refinement_increases_both_regions_by_eight():
    a=geometry_metrics(geometry_profile('paper','coarse'))
    b=geometry_metrics(geometry_profile('paper','medium'))
    c=geometry_metrics(geometry_profile('paper','fine'))
    for phase in ('cells_fluid','cells_copper','cells_total'):
        assert b[phase]==8*a[phase]
        assert c[phase]==8*b[phase]


def test_water_energy_and_channel_pressure_baseline():
    g=geometry_profile()
    c=experimental_condition(890,100,g)
    assert 0.5<c['flow_L_min']<0.8
    assert 7<c['adiabatic_outlet_rise_K']<15
    assert math.isclose(c['massflow_kg_s']*4180*c['adiabatic_outlet_rise_K'],c['power_W'],rel_tol=1e-12)
    assert laminar_core_pressure_pa(1000)>laminar_core_pressure_pa(500)>0
    # CFD channel-only must not be presented as measured plenum-to-plenum.
    assert 0<laminar_core_pressure_pa(890)<86_000


def test_mesh_case_has_no_fabricated_output(tmp_path):
    g=geometry_profile()
    folder=tmp_path/'case'
    write_geometry_case(g,folder)
    write_fv_dictionaries(folder,g,890,100)
    write_thermophysics(folder)
    assert (folder/'system/blockMeshDict').exists()
    assert not (folder/'0/water/T').exists() # until actual split patch names are known
    assert not (folder/'log.cht').exists()


def test_comparison_rejects_missing_openfoam_results(tmp_path):
    with pytest.raises(FileNotFoundError,match='NOT RUN'):
        compare_against_paper(tmp_path/'cfd_export.csv',Path(__file__).resolve().parents[1])


def test_foam_boundary_writer_requires_actual_interface(tmp_path):
    from openchipcool.cfd.foam_setup import write_fields_after_split
    w=tmp_path/'constant/water/polyMesh'
    c=tmp_path/'constant/copper/polyMesh'
    w.mkdir(parents=True);c.mkdir(parents=True)
    w.joinpath('boundary').write_text('''FoamFile
{
 version 2.0;
}
4
(
inlet
{
 type patch;
 nFaces 1;
}
outlet
{
 type patch;
 nFaces 1;
}
cover
{
 type wall;
 nFaces 1;
}
water_to_copper
{
 type mappedWall;
 nFaces 1;
}
);
''')
    c.joinpath('boundary').write_text('''3
(
heater
{
 type wall;
 nFaces 1;
}
copperEnds
{
 type wall;
 nFaces 1;
}
copper_to_water
{
 type mappedWall;
 nFaces 1;
}
);
''')
    created=write_fields_after_split(tmp_path,2.5,288.15,100.)
    assert len(created)==6
    assert 'fixedGradient' in (tmp_path/'0/copper/T').read_text()
    assert 'compressible::turbulentTemperatureCoupledBaffleMixed' in (tmp_path/'0/water/T').read_text()
    assert 'fixedValue' in (tmp_path/'0/water/U').read_text()


def test_blockface_audit(tmp_path):
    from openchipcool.cfd.audit import audit_mesh_blocks
    for name in ('paper','design36'):
        data=audit_mesh_blocks(geometry_profile(name,'coarse'))
        assert data['copper_water_interface_block_faces']==3*geometry_profile(name,'coarse').n_channels
        assert data['patches']['inlet']==geometry_profile(name,'coarse').n_channels
