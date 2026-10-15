import numpy as np
import pytest
from openchipcool.config import Spec
from openchipcool.fem import build_fem,run_fem
from openchipcool.scenarios import make_times,training_power

@pytest.fixture(scope="module")
def fem():
    return build_fem(nx=10)

def test_mesh_material_and_source_conservation(fem):
    s=Spec()
    assert fem.coords.shape[1]==3
    assert fem.tets.shape[1]==4
    assert np.all(fem.C>0)
    assert np.allclose(fem.Q.sum(axis=0),np.ones(4),atol=1e-12)
    assert np.isclose(np.sum(fem.area),s.side**2,atol=1e-14)
    assert np.isclose(sum(fem.volume_by_material.values()),s.side**2*s.height,rtol=1e-10)
    assert np.allclose((fem.K-fem.K.T).data,0,atol=1e-10)
    assert np.max(np.abs(fem.K@np.ones(len(fem.C))))<1e-8


def test_steady_energy_balance(fem):
    power=[200,180,220,160]
    T=fem.steady(power,inlet=30,h=26000)
    assert np.all(np.isfinite(T))
    assert np.all(fem.average_die_temps(T)>30)
    assert np.isclose(fem.boundary_heat(T,30,26000),sum(power),rtol=1e-9)
    assert np.max(fem.peak_die_temps(T)) >= np.max(fem.average_die_temps(T))


def test_transient_implicit_energy_balance(fem):
    times=[0.,0.1]
    power=np.array([[0.,0.,0.,0.],[150.,180.,140.,200.]])
    result=run_fem(fem,times,power,inlet=30,h=26000)
    u=result["last_field"]
    dE=np.dot(fem.C,u-30)/.1
    q=fem.boundary_heat(u,30,26000)
    assert np.isclose(dE+q,np.sum(power[-1]),rtol=1e-9)


def test_transient_shape(fem):
    t=make_times(1.5,.15)
    r=run_fem(fem,t,training_power(t))
    assert r["means"].shape==(len(t),4)
    assert r["peaks"].shape==(len(t),4)
    assert np.allclose(r["means"][0],30.)
