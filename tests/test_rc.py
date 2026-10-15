import numpy as np
from openchipcool.config import Spec,h_from_flow
from openchipcool.fem import build_fem,run_fem
from openchipcool.rc import from_geometry,fit_to_fem,run_rc_fixed_inlet,simulate_loop
from openchipcool.scenarios import make_times,training_power,holdout_power,dashboard_power
from openchipcool.modelica import render_modelica


def test_positive_components_and_flow_law():
    r=from_geometry()
    assert min(r.C_die,r.C_plate,r.G_vertical,r.G_lateral,r.G_conv_ref)>0
    assert h_from_flow(0,Spec())==0
    assert h_from_flow(.2,Spec())<h_from_flow(4,Spec())<h_from_flow(8,Spec())


def test_holdout_matches_fem_means_on_reference_mesh():
    f=build_fem(nx=10)
    t=make_times()
    train=run_fem(f,t,training_power(t))
    fitted=fit_to_fem(t,training_power(t),train["means"])
    valid=run_fem(f,t,holdout_power(t))
    pred=run_rc_fixed_inlet(t,holdout_power(t),fitted)
    rmse=np.sqrt(np.mean((pred[:,:4]-valid["means"])**2))
    assert fitted.calibrated
    assert rmse<3.0  # coarse-mesh regression bound, not experimental accuracy claim
    assert np.max(valid["peaks"])>np.max(valid["means"])


def test_pump_fault_increases_die_mean():
    s=Spec();r=from_geometry()
    t=np.arange(0,60.01,.5)
    p=dashboard_power(t,800,"step")
    normal=simulate_loop(t,p,np.full(len(t),4.),r,s)
    fault=simulate_loop(t,p,np.where(t<25,4.,.2),r,s)
    assert fault["die"][-1].max()>normal["die"][-1].max()+5
    assert np.all(np.isfinite(fault["states"]))


def test_modelica_export_parameters():
    src=render_modelica(from_geometry())
    for marker in ["model OpenChipCoolLoop", "der(Tdie[i])", "der(Tcool)", "pumpFailureTime", "GconvRef", "nusseltCurrent", "ReCurrent"]:
        assert marker in src
