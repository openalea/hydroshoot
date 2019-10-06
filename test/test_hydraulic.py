from numpy import arange
from pytest import fixture

from openalea.mtg import traversal

from hydroshoot import hydraulic, architecture
from non_regression_data import potted_syrah


@fixture()
def simple_shoot():
    g = potted_syrah()
    vid_base = architecture.mtg_base(g, vtx_label='inT')
    g.node(g.root).vid_base = vid_base
    return g


def test_conductivity_max_increases_as_segment_diameter_increases():
    cond_max = [hydraulic.conductivity_max(diameter, a=2.8, b=0.1, min_kmax=0.) for diameter in arange(0, 0.2, 0.005)]
    assert all(x < y for x, y in zip(cond_max, cond_max[1:]))


def test_conductivity_max_is_higher_or_equal_to_its_minimal_value():
    min_kmax = 0.02
    assert hydraulic.conductivity_max(diameter=0., a=2.8, b=0.1, min_kmax=min_kmax) >= min_kmax


def test_conductivity_max_changes_following_parameters_values():
    cond_max_a = [hydraulic.conductivity_max(diameter=0.01, a=a, b=0.1, min_kmax=0.) for a in arange(2.5, 2.8, 0.1)]
    cond_max_b = [hydraulic.conductivity_max(diameter=0.01, a=2.8, b=b, min_kmax=0.) for b in arange(2., 5., 0.5)]
    assert all(x < b for x, y in zip(cond_max_a, cond_max_a[1:]))
    assert all(x < b for x, y in zip(cond_max_b, cond_max_b[1:]))


def test_cavitation_factor_raises_error_if_undefined_model():
    try:
        hydraulic.cavitation_factor(psi=0.0, model='some_model', fifty_cent=-0.51, sig_slope=3)
    except ValueError as err:
        assert err.args[0] == "The 'model' argument must be one of the following ('misson','tuzet', 'linear')."


def test_cavitation_factor_increases_as_segment_water_potential_decreases():
    for model in ('misson', 'tuzet', 'linear'):
        cavitation = [hydraulic.cavitation_factor(psi, model, -0.51, 3) for psi in arange(0, -3, -0.1)]
        assert all(x >= y for x, y in zip(cavitation, cavitation[1:]))


def test_def_param_soil_returns_the_right_soil_property_values():
    ref_values = {'Sand': (0.045, 0.430, 0.145, 2.68, 712.8),
                  'Loamy_Sand': (0.057, 0.410, 0.124, 2.28, 350.2),
                  'Sandy_Loam': (0.065, 0.410, 0.075, 1.89, 106.1),
                  'Loam': (0.078, 0.430, 0.036, 1.56, 24.96),
                  'Silt': (0.034, 0.460, 0.016, 1.37, 6.00),
                  'Silty_Loam': (0.067, 0.450, 0.020, 1.41, 10.80),
                  'Sandy_Clay_Loam': (0.100, 0.390, 0.059, 1.48, 31.44),
                  'Clay_Loam': (0.095, 0.410, 0.019, 1.31, 6.24),
                  'Silty_Clay_Loam': (0.089, 0.430, 0.010, 1.23, 1.68),
                  'Sandy_Clay': (0.100, 0.380, 0.027, 1.23, 2.88),
                  'Silty_Clay': (0.070, 0.360, 0.005, 1.09, 0.48),
                  'Clay': (0.068, 0.380, 0.008, 1.09, 4.80)}
    module_values = hydraulic.def_param_soil()
    for soil_class, soil_properties in ref_values.items():
        assert all(x == y for x, y in zip(soil_properties, module_values[soil_class]))


def test_def_param_soil_returns_costum_soil_properties_when_provided():
    costum_soil = (0.02, 0.3, 0.03, 1.5, 25)
    assert hydraulic.def_param_soil(costum_soil)['Costum'] == costum_soil


def test_k_soil_soil_decreases_as_water_potential_decreases():
    for soil_class in hydraulic.def_param_soil().keys():
        soil_conductivity = [hydraulic.k_soil_soil(psi, soil_class) for psi in arange(0, -3, -0.1)]
        assert all(x >= y for x, y in zip(soil_conductivity, soil_conductivity[1:]))


def test_k_soil_soil_maximum_value_is_greater_for_sand_than_clay():
    soil_conductivity = [hydraulic.k_soil_soil(0., soil_class) for soil_class in ('Clay', 'Sand')]
    assert all(x <= y for x, y in zip(soil_conductivity, soil_conductivity[1:]))


def k_soil_root_increases_as_soil_conductivity_increases():
    soil_conductivity = [0.48, 1.68, 2.88, 4.8, 6.0, 6.24, 10.8, 24.96, 31.44, 106.1, 350.2, 712.8]
    soil_root_cond = [hydraulic.k_soil_root(k_soil, dist_roots=0.013, rad_roots=.0001) for k_soil in soil_conductivity]
    assert all(x <= y for x, y in zip(soil_root_cond, soil_root_cond[1:]))


def k_soil_root_increases_as_the_ratio_of_average_root_distance_to_radius_decreases():
    soil_root_cond = [hydraulic.k_soil_root(k_soil=6., dist_roots=0.013, rad_roots=0.013 / x)
                      for x in range(200, 100, -10)]
    assert all(x <= y for x, y in zip(soil_root_cond, soil_root_cond[1:]))


def test_soil_water_potential_decreases_as_water_withdrawal_increases():
    for soil_class in hydraulic.def_param_soil().keys():
        psi_soil = [hydraulic.soil_water_potential(psi_soil_init=0., water_withdrawal=w, soil_class=soil_class,
                                                   soil_total_volume=1, psi_min=-3.)
                    for w in arange(0, 3, 0.1)]
        assert all(x >= y for x, y in zip(psi_soil, psi_soil[1:]))


def test_soil_water_potential_drops_faster_for_small_soil_reservoirs_than_bigger_ones():
    for soil_class in hydraulic.def_param_soil().keys():
        psi_soil_small = hydraulic.soil_water_potential(psi_soil_init=0., water_withdrawal=1., soil_class=soil_class,
                                                        soil_total_volume=1., psi_min=-3.)
        psi_soil_big = hydraulic.soil_water_potential(psi_soil_init=0., water_withdrawal=1., soil_class=soil_class,
                                                      soil_total_volume=2., psi_min=-3.)
        assert psi_soil_small < psi_soil_big


def test_hydraulic_prop_attributes_the_right_properties_to_the_different_organs(simple_shoot):
    vid_base = simple_shoot.node(simple_shoot.root).vid_base
    for vtx_id in traversal.post_order2(simple_shoot, vid_base):
        n = simple_shoot.node(vtx_id)
        if n.label.startswith('LI'):
            n.E = 0.
            n.An = 0.

    hydraulic.hydraulic_prop(simple_shoot, mass_conv=18.01528, length_conv=1.e-2, a=2.6, b=2.0, min_kmax=0.)

    for vtx_id in traversal.post_order2(simple_shoot, vid_base):
        n = simple_shoot.node(vtx_id)
        assert hasattr(n, 'Flux')
        assert hasattr(n, 'FluxC')
        if n.label.startswith(('in', 'cx', 'Pet')):
            assert hasattr(n, 'Kmax')

# def test_xylem_water_potential(simple_shoot):
#     vid_base = simple_shoot.node(simple_shoot.root).vid_base
#     for vtx_id in traversal.post_order2(simple_shoot, vid_base):
#         n = simple_shoot.node(vtx_id)
#         n.psi_head = 0.
#         if n.label.startswith('LI'):
#             n.E = 0.01
#             n.An = 0.
#
#     simple_shoot.node(vid_base).psi_head = 0.
#     hydraulic.hydraulic_prop(simple_shoot, mass_conv=18.01528, length_conv=1.e-2, a=2.6, b=2.0, min_kmax=0.)
#
#
#     hydraulic.xylem_water_potential(simple_shoot, psi_soil=-0.8, model='tuzet', psi_min=-3.0, psi_error_crit=0.001, max_iter=100,
#                           length_conv=1.E-2, fifty_cent=-0.51, sig_slope=0.1, dist_roots=0.013, rad_roots=.0001,
#                           negligible_shoot_resistance=False, start_vid=None, stop_vid=None, psi_step=0.5)
