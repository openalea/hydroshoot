from numpy import arange
from pytest import approx, fixture

from hydroshoot import exchange


@fixture()
def photosynthesis_parameters():
    return exchange.par_photo_default()


def test_leaf_na_is_as_expected():
    expected_result = 2.007
    obtained_result = exchange.leaf_Na(age_gdd=1000., ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635)

    assert obtained_result == approx(expected_result, abs=1.e-3)


def test_leaf_na_reduces_as_age_gdd_increases():
    assert exchange.leaf_Na(
        age_gdd=0., ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635) > exchange.leaf_Na(
        age_gdd=1000., ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635)


def test_leaf_na_reduces_as_cumulative_ppfd_reduces():
    assert exchange.leaf_Na(
        age_gdd=1000., ppfd_10=0.0, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635) < exchange.leaf_Na(
        age_gdd=1000., ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635)


def test_leaf_na_is_null_when_cumulative_ppfd_and_relationship_intercept_are_null():
    assert exchange.leaf_Na(age_gdd=1000., ppfd_10=0.0, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=0.0) == 0.0


def test_leaf_na_is_null_when_age_gdd_and_relationship_intercept_are_null():
    assert exchange.leaf_Na(age_gdd=0.0, ppfd_10=38.64, a_n=-0.0008, b_n=0.0, a_m=6.471, b_m=56.635) == 0.0


def test_leaf_na_is_greater_or_equal_to_zero():
    assert exchange.leaf_Na(age_gdd=10000.0, ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635) == 0.0


def test_arrhenius_1_increases_as_temperature_increases(photosynthesis_parameters):
    param_names = ['Tx', 'Kc', 'Ko']
    for param_name in param_names:
        prev_value = exchange.arrhenius_1(param_name, 0, photosynthesis_parameters)
        for leaf_temperature in range(1, 50):
            actual_value = exchange.arrhenius_1(param_name, leaf_temperature, photosynthesis_parameters)
            assert actual_value > prev_value
            prev_value = actual_value


def test_arrhenius_2_is_maximum_for_vcmax_at_39_degrees_celsius(photosynthesis_parameters):
    value_at_low_temperature = exchange.arrhenius_2('Vcmax', 0.0, photosynthesis_parameters)
    value_at_optimal_temperature = exchange.arrhenius_2('Vcmax', 39.0, photosynthesis_parameters)
    value_at_high_temperature = exchange.arrhenius_2('Vcmax', 50.0, photosynthesis_parameters)

    assert all([value_at_optimal_temperature > val for val in (value_at_low_temperature, value_at_high_temperature)])


def test_arrhenius_2_is_maximum_for_jmax_at_37_degrees_celsius(photosynthesis_parameters):
    value_at_low_temperature = exchange.arrhenius_2('Jmax', 0.0, photosynthesis_parameters)
    value_at_optimal_temperature = exchange.arrhenius_2('Jmax', 37.0, photosynthesis_parameters)
    value_at_high_temperature = exchange.arrhenius_2('Jmax', 50.0, photosynthesis_parameters)

    assert all([value_at_optimal_temperature > val for val in (value_at_low_temperature, value_at_high_temperature)])


def test_dHd_sensibility_decreases_as_leaf_water_potential_decreases():
    prev_value = exchange.dHd_sensibility(psi=0.0,
                                          temp=25, dhd_max=200.,
                                          dhd_inhib_beg=195., dHd_inhib_max=190.,
                                          psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                          temp_inhib_beg=35, temp_inhib_max=40)
    for psi in arange(0, -3, -0.1):
        act_value = exchange.dHd_sensibility(psi,
                                             temp=25, dhd_max=200.,
                                             dhd_inhib_beg=195., dHd_inhib_max=190.,
                                             psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                             temp_inhib_beg=35, temp_inhib_max=40)
        assert act_value <= prev_value
        prev_value = act_value


def test_dHd_sensibility_is_not_affected_by_leaf_temperature_increases_when_leaf_water_potential_is_optimal():
    value_at_optimal_psi = exchange.dHd_sensibility(psi=0.0,
                                                    temp=25, dhd_max=200.,
                                                    dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                    psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                    temp_inhib_beg=35, temp_inhib_max=40)
    for temperature in range(0, 50):
        assert value_at_optimal_psi == exchange.dHd_sensibility(psi=0.0,
                                                                temp=temperature, dhd_max=200.,
                                                                dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                                psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                                temp_inhib_beg=35, temp_inhib_max=40)


def test_dHd_sensibility_is_affected_by_leaf_temperature_when_leaf_water_potential_is_below_its_minimum_threshold():
    value_at_optimal_psi = exchange.dHd_sensibility(psi=0.0,
                                                    temp=25, dhd_max=200.,
                                                    dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                    psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                    temp_inhib_beg=35, temp_inhib_max=40)
    for temperature in range(0, 50):
        assert value_at_optimal_psi > exchange.dHd_sensibility(psi=-2.0,
                                                               temp=temperature, dhd_max=200.,
                                                               dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                               psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                               temp_inhib_beg=35, temp_inhib_max=40)


def test_dHd_sensibility_is_affected_by_leaf_temperature_only_within_restricted_bounds():
    low_temperature_bound = 35
    high_temperature_bound = 40
    value_at_lowest_bound = exchange.dHd_sensibility(psi=-2.0,
                                                     temp=low_temperature_bound, dhd_max=200.,
                                                     dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                     psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                     temp_inhib_beg=35, temp_inhib_max=40)
    value_at_highest_bound = exchange.dHd_sensibility(psi=-2.0,
                                                      temp=high_temperature_bound, dhd_max=200.,
                                                      dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                      psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                      temp_inhib_beg=35, temp_inhib_max=40)
    for temperature in range(0, 50):
        if temperature < low_temperature_bound:
            assert value_at_lowest_bound == exchange.dHd_sensibility(psi=-2.0,
                                                                     temp=temperature, dhd_max=200.,
                                                                     dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                                     psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                                     temp_inhib_beg=35, temp_inhib_max=40)
        elif temperature < high_temperature_bound:
            assert value_at_lowest_bound >= exchange.dHd_sensibility(psi=-2.0,
                                                                     temp=temperature, dhd_max=200.,
                                                                     dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                                     psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                                     temp_inhib_beg=35, temp_inhib_max=40)
        else:
            assert value_at_highest_bound == exchange.dHd_sensibility(psi=-2.0,
                                                                      temp=temperature, dhd_max=200.,
                                                                      dhd_inhib_beg=195., dHd_inhib_max=190.,
                                                                      psi_inhib_beg=-.75, psi_inhib_max=-2.,
                                                                      temp_inhib_beg=35, temp_inhib_max=40)


def test_compute_an_2par_increases_only_electron_transport_as_ppfd_increases(photosynthesis_parameters):
    vcmax = []
    j_frac = []
    tpu_triple = []
    rd = []

    for ppfd in range(0, 2000, 10):
        res = exchange.compute_an_2par(photosynthesis_parameters, ppfd, leaf_temp=25.0)
        vcmax.append(res[0])
        j_frac.append(res[2])
        tpu_triple.append(res[4])
        rd.append(res[6])

    # vcmax is independent from ppfd
    assert all(x == y for x, y in zip(vcmax, vcmax[1:]))
    # j increases with ppfd
    assert all(x < y for x, y in zip(j_frac, j_frac[1:]))
    # tpu is independent from ppfd
    assert all(x == y for x, y in zip(tpu_triple, tpu_triple[1:]))
    # rd is independent from ppfd
    assert all(x == y for x, y in zip(rd, rd[1:]))


def test_compute_an_2par_affects_all_photosynthetic_parameters_by_temperature(photosynthesis_parameters):
    vcmax = []
    j_frac = []
    tpu_triple = []
    rd = []

    for temperature in range(-10, 45):
        res = exchange.compute_an_2par(photosynthesis_parameters, ppfd=1800., leaf_temp=temperature)
        vcmax.append(res[0])
        j_frac.append(res[2])
        tpu_triple.append(res[4])
        rd.append(res[6])

    # vcmax is dependent on temperature
    assert all(x != y for x, y in zip(vcmax, vcmax[1:]))
    # j is dependent on temperature
    assert all(x != y for x, y in zip(j_frac, j_frac[1:]))
    # tpu is dependent on temperature
    assert all(x != y for x, y in zip(tpu_triple, tpu_triple[1:]))
    # rd is dependent on temperature
    assert all(x != y for x, y in zip(rd, rd[1:]))


def test_fvpd_3_raises_error_for_unrecognized_model():
    try:
        exchange.fvpd_3('any_model', vpd=3., psi=-0., psi_crit=-0.37, m0=5.278, steepness_tuzet=1.85, d0_leuning=30.)
    except ValueError as err:
        assert err.args[0] == "The 'model' argument must be one of the following ('misson','tuzet', 'linear' or 'vpd')."


def test_fvpd_3_increases_stress_factor_as_leaf_water_potential_decreases():
    for model in ('misson', 'tuzet', 'linear', 'vpd'):
        reduction_factor = []
        for psi in arange(0, -3, -0.1):
            reduction_factor.append(
                exchange.fvpd_3(model, vpd=3., psi=psi, psi_crit=-0.37, m0=5.278, steepness_tuzet=1.85, d0_leuning=30.))
        assert all(x >= y for x, y in zip(reduction_factor, reduction_factor[1:]))


def test_mesophyll_conductance_changes_with_temperature():
    gm = [exchange.mesophyll_conductance(t, gm25=0.1025, activation_energy=49600., deactivation_energy=437400.,
                                         entropy=1400.) for t in arange(-10, 46)]
    assert all(x != y for x, y in zip(gm, gm[1:]))


def test_mesophyll_conductance_is_positive():
    gm = [exchange.mesophyll_conductance(t, gm25=0.1025, activation_energy=49600., deactivation_energy=437400.,
                                         entropy=1400.) for t in arange(-10, 46)]
    assert all(x >= 0 for x in gm)


def test_mesophyll_conductance_is_maximum_at_36_degrees_celsius():
    gm_max = exchange.mesophyll_conductance(36, gm25=0.1025, activation_energy=49600., deactivation_energy=437400.,
                                            entropy=1400.)
    gm = [exchange.mesophyll_conductance(t, gm25=0.1025, activation_energy=49600., deactivation_energy=437400.,
                                         entropy=1400.) for t in arange(-10, 46)]
    assert all(x <= gm_max for x in gm)


def test_boundary_layer_conductance_decreases_as_leaf_length_increases():
    gb = [exchange.boundary_layer_conductance(l, wind_speed=2.0, atm_pressure=101.3, air_temp=25.,
                                              ideal_gas_cst=exchange.R) for l in arange(0.001, 0.25, 0.02)]
    assert all(x >= y for x, y in zip(gb, gb[1:]))


def test_boundary_layer_conductance_increases_as_wind_speed_increases():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=w, atm_pressure=101.3, air_temp=25.,
                                              ideal_gas_cst=exchange.R) for w in arange(0.0, 5.0, 0.5)]
    assert all(x <= y for x, y in zip(gb, gb[1:]))


def test_boundary_layer_conductance_is_weakly_dependent_on_atmospheric_pressure():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=2., atm_pressure=p, air_temp=25.,
                                              ideal_gas_cst=exchange.R) for p in arange(90.3, 102.3)]
    assert all(x == approx(y, abs=1.e-3) for x, y in zip(gb, gb[1:]))


def test_boundary_layer_conductance_increases_as_air_temperature_increases():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=2.0, atm_pressure=101.3, air_temp=t,
                                              ideal_gas_cst=exchange.R) for t in range(-10, 46)]
    assert all(x <= y for x, y in zip(gb, gb[1:]))
