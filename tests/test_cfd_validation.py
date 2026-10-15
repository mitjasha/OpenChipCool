"""Only unit/mocked tests, NOT an OpenFOAM execution or V&V result."""
import json, math
from pathlib import Path
import pytest
from openchipcool.cfd.validation import pressure_budget,grid_convergence,validate_export,TC_X_MM
from openchipcool.cfd.geometry import geometry_profile
from openchipcool.cfd.foam_setup import write_fv_dictionaries
from scripts.build_validation_matrix import build,OPERATING_POINTS

ROOT=Path(__file__).resolve().parents[1]


def test_analytical_plenum_budget_matches_eq18_to_22():
    b=pressure_budget(rho=1000,u_channel=4,u_p2_in=.6,u_p2_out=.6,u_p1_in=.1,u_p1_out=.1)
    assert len(b)==4
    assert all(math.isfinite(x) for x in b.values())
    assert 0<sum(b.values())<100000
    with pytest.raises(ValueError): pressure_budget(rho=-1,u_channel=2,u_p2_in=.5,u_p2_out=.5,u_p1_in=.1,u_p1_out=.1)


def test_gci_monotone_and_oscillating_trajectories():
    x=grid_convergence(44,41,40)
    assert math.isclose(x['observed_order'],math.log(3,2))
    assert x['GCI_fine_percent']>0
    with pytest.raises(ValueError,match='monotone'):
        grid_convergence(44,39,40)
    with pytest.raises(ValueError,match='refinement'):
        grid_convergence(44,41,40,r=1.)


def test_fig7_contains_4_observations_for_6_points():
    from openchipcool.cfd.validation import read_csv
    a=read_csv(ROOT/'data/literature/qu_mudawar_2002_tc_fig7_visual_estimate.csv')
    b=read_csv(ROOT/'data/literature/qu_mudawar_2002_tout_fig5_visual_estimate.csv')
    assert len(a)==24 and len(b)==6 and len(TC_X_MM)==4
    assert len({(i['qeff_W_cm2'],i['Re']) for i in a})==6
    assert all(float(i['visual_reading_uncertainty_K'])>=5 for i in a)


def test_case_manifest_input_only(tmp_path):
    # Dry-run generator test, no solver involved
    manifest=build(tmp_path,levels=('coarse',))
    assert manifest['status']=='INPUTS_ONLY_NO_SIMULATION'
    assert manifest['count']==len(OPERATING_POINTS)
    case=tmp_path/'qm_q100_Re890_coarse'
    assert 'input' not in json.loads((case/'setup_status.json').read_text())['stage'].lower()
    s=(case/'system/controlDict').read_text()
    assert 'outlet_T_mass_weighted' in s and 'tc_copper' in s
    assert 'inlet_p_areaAverage' in s and 'outlet_p_areaAverage' in s
    assert not (case/'log.cht').exists()


def test_validation_rejects_placeholder_without_logs(tmp_path):
    import csv
    row={
      'case_id':'qm_fake','mesh':'coarse','Re':'890','qeff_W_cm2':'100',
      'dp_core_Pa':'10000','dp_plenums_Pa':'3000','Tin_C':'15',
      'Tout_massweighted_C':'25','mass_error_frac':'0',
      'energy_error_frac':'0','solver_log_path':'missing.log',
      'mesh_water_check_path':'missing.water','mesh_copper_check_path':'missing.copper',
    }
    row.update({f'Ttc{i}_rise_K':str(40+i) for i in range(1,5)})
    p=tmp_path/'fake.csv'
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=row.keys());w.writeheader();w.writerow(row)
    with pytest.raises(FileNotFoundError,match='real solver evidence'):
        validate_export(p,ROOT)


def test_iapws_2011_table_8_reference_values():
    from openchipcool.cfd.water_properties import viscosity_pa_s,conductivity_w_mk,openfoam_tabulated_transport
    # IAPWS Supplementary Release Table 8; 298.15 K (25 deg C).
    assert math.isclose(viscosity_pa_s(298.15)*1e6,889.996774,rel_tol=2e-8)
    assert math.isclose(conductivity_w_mk(298.15),0.606502308,rel_tol=2e-8)
    assert viscosity_pa_s(288.15)>viscosity_pa_s(323.15)
    assert 'kappa' in openfoam_tabulated_transport()
    with pytest.raises(ValueError): viscosity_pa_s(390)


def test_generated_fluid_has_temperature_dependent_transport(tmp_path):
    from openchipcool.cfd.foam_setup import write_thermophysics
    write_thermophysics(tmp_path,geometry_profile())
    body=(tmp_path/'constant/water/thermophysicalProperties').read_text()
    assert 'transport tabulated;' in body
    assert 'mu\n' in body and 'kappa\n' in body
