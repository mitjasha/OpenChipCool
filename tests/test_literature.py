import numpy as np
from openchipcool.config import Spec,h_from_flow
from openchipcool.literature import ChannelDesign,reynolds,reference_nusselt,equivalent_htc
from openchipcool.fem import build_aligned_fem
from openchipcool.cli import richardson_gci


def test_geometry_exact_across_refinements():
    s=Spec()
    for r in (1,2):
        f=build_aligned_fem(s,refinement=r)
        assert np.allclose(f.source_volume,s.chip_side**2*s.silicon_thickness,rtol=1e-12)
        assert np.isclose(np.sum(f.area),s.side**2,rtol=1e-12)
        assert np.allclose(np.sum(f.Q,axis=0),1,atol=1e-12)
        assert np.allclose(f.observation.sum(axis=1),1,atol=1e-12)
        assert np.max(abs(f.boundary_heat(f.steady(h=s.h_ref),s.coolant_initial,s.h_ref)-800)) <1e-5


def test_literature_reference_and_applicability():
    c=ChannelDesign();b,side=reference_nusselt()
    assert 4<b<9 and 7<side<13
    assert 139 < reynolds(4,c) < 1672
    assert h_from_flow(0,Spec())==0
    assert 15000 < equivalent_htc(4,c) < 150000
    assert equivalent_htc(.2,c) < equivalent_htc(4,c)


def test_gci_known_second_order_and_oscillatory_flag():
    refs=np.array([1,2,3])
    data=10+2/refs**2
    g=richardson_gci(data,tuple(refs))
    assert np.isclose(g['observed_order'],2,rtol=1e-7)
    assert np.isclose(g['extrapolated_value'],10,rtol=1e-10)
    assert richardson_gci([1,3,2])['status']=='inconclusive'
