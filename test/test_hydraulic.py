from numpy import arange
from openalea.mtg import traversal

from openalea.hydroshoot import hydraulic, architecture
from non_regression_data import potted_syrah


def test_conductivity_max_increases_as_segment_diameter_increases():
    cond_max = [hydraulic.conductivity_max(diameter, a=2.8, b=0.1, min_kmax=0.) for diameter in arange(0, 0.2, 0.005)]
    assert all(x < y for x, y in zip(cond_max, cond_max[1:]))


def test_conductivity_max_is_higher_or_equal_to_its_minimal_value():
    min_kmax = 0.02
    assert hydraulic.conductivity_max(diameter=0., a=2.8, b=0.1, min_kmax=min_kmax) >= min_kmax


def test_conductivity_max_changes_following_parameters_values():
    cond_max_a = [hydraulic.conductivity_max(diameter=0.01, a=a, b=0.1, min_kmax=0.) for a in arange(2.5, 2.8, 0.1)]
    cond_max_b = [hydraulic.conductivity_max(diameter=0.01, a=2.8, b=b, min_kmax=0.) for b in arange(2., 5., 0.5)]
    assert all(x < y for x, y in zip(cond_max_a, cond_max_a[1:]))
    assert all(x > y for x, y in zip(cond_max_b, cond_max_b[1:]))


def test_cavitation_factor_raises_error_if_undefined_model():
    try:
        hydraulic.cavitation_factor(psi=0.0, model='some_model', fifty_cent=-0.51, sig_slope=3)
    except ValueError as err:
        assert err.args[0] == "The 'model' argument must be one of the following ('misson','tuzet', 'linear')."


def test_cavitation_factor_increases_as_segment_water_potential_decreases():
    for model in ('misson', 'tuzet', 'linear'):
        cavitation = [hydraulic.cavitation_factor(psi, model, -0.51, 3) for psi in arange(0, -3, -0.1)]
        assert all(x >= y for x, y in zip(cavitation, cavitation[1:]))


def test_hydraulic_prop_attributes_the_right_properties_to_the_different_organs():
    simple_shoot = potted_syrah()
    simple_shoot.node(simple_shoot.root).vid_base = architecture.get_mtg_base(simple_shoot, vtx_label='inT')

    vid_base = simple_shoot.node(simple_shoot.root).vid_base
    for vtx_id in traversal.post_order2(simple_shoot, vid_base):
        n = simple_shoot.node(vtx_id)
        if n.label.startswith('LI'):
            n.E = 0.
            n.An = 0.

    hydraulic.hydraulic_prop(simple_shoot, length_conv=1.e-2, a=2.6, b=2.0, min_kmax=0.)

    for vtx_id in traversal.post_order2(simple_shoot, vid_base):
        n = simple_shoot.node(vtx_id)
        assert hasattr(n, 'Flux')
        assert hasattr(n, 'FluxC')
        if n.label.startswith(('in', 'cx', 'Pet')):
            assert hasattr(n, 'Kmax')
