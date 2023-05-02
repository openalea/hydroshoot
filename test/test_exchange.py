from datetime import datetime

from numpy import arange, linspace, testing
from pandas import Series

from hydroshoot import exchange, utilities, constants


def setup_default_photosynthesis_params():
    """
    Generates a dictionary containing default **25 degrees** paramter values of the Farquhar's model for Vitis vinifera cv. Syrah.

   :Parameters:
   - **Vcm25**: Maximum RuBP-saturated rate of carboxylation [umol m-2 s-1]
   - **Jm25**: Maximum of electron transport [umol m-2 s-1]
   - **cRd**: the coefficient of mitochondrial respiration to `Vcm25`
   - **TPU25**: the rate of triose phosphate transport [umol m-2 s-1]
   - **Kc25**: Michaelis-Menten constant for the carboxylase [umol mol-1]
   - **Ko25**: Michaelis-Menten constant for the oxygenase [mmol mol-1]
   - **Tx25**: CO2 compensation point in the absence of mitochondrial respiration [umol mol-1]
   - **c**: empirical parameter defining the temperature response curves of each of Kc, Ko, Vcm, Jm, TPU, Rd and Tx
   - **deltaHa**: Activation energy of the Arrhenius functions [kJ molCO2-1]
   - **ds**: float, enthalpy of activation [KJ mol-1]
   - **dHd**: float, enthalpy of deactivation [KJ mol-1]
   """

    res = dict(
        Vcm25=89.0,
        Jm25=143.0,
        cRd=0.008,
        TPU25=10.0,
        alpha=0.24,
        Kc25=404.9,
        Ko25=278.4,
        Tx25=42.75,
        ds=0.635,
        dHd=200,
        RespT_Kc={'model': 'Ahrenius', 'c': 38.05, 'deltaHa': 79.43},
        RespT_Ko={'model': 'Ahrenius', 'c': 20.30, 'deltaHa': 36.38},
        RespT_Vcm={'model': 'Ahrenius', 'c': 26.35, 'deltaHa': 65.33},
        RespT_Jm={'model': 'Ahrenius', 'c': 17.57, 'deltaHa': 43.54},
        RespT_TPU={'model': 'Ahrenius', 'c': 21.46, 'deltaHa': 53.1},
        RespT_Rd={'model': 'Ahrenius', 'c': 18.72, 'deltaHa': 46.39},
        RespT_Tx={'model': 'Ahrenius', 'c': 19.02, 'deltaHa': 37.83}
    )
    res.update({'Rd': res['cRd'] * res['Vcm25']})
    return res


def setup_default_stomatal_conductance_params():
    return dict(model='misson', g0=0.019, m0=5.278, psi0=-0.1, D0=30., n=1.85)


def setup_leaf_local_weather():
    return Series({'time': datetime(2012, 8, 1, 11),
                   'Tac': 26.84,
                   'hs': 43.77,
                   'Rg': 300.28742,
                   'u': 2.17,
                   'Ca': 400,
                   'Pa': 101.3,
                   'PPFD': 663.035})


def test_leaf_na_is_as_expected():
    expected_result = 2.007
    obtained_result = exchange.leaf_Na(age_gdd=1000., ppfd_10=38.64, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635)

    testing.assert_almost_equal(obtained_result, expected_result, decimal=3)


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


def test_arrhenius_1_increases_as_temperature_increases():
    param_names = ['Tx', 'Kc', 'Ko']
    for param_name in param_names:
        prev_value = exchange.arrhenius_1(param_name, 0, setup_default_photosynthesis_params())
        for leaf_temperature in range(1, 50):
            actual_value = exchange.arrhenius_1(param_name, leaf_temperature,
                                                setup_default_photosynthesis_params())
            assert actual_value > prev_value
            prev_value = actual_value


def test_arrhenius_2_is_maximum_for_vcmax_at_39_degrees_celsius():
    value_at_low_temperature = exchange.arrhenius_2('Vcmax', 0.0, setup_default_photosynthesis_params())
    value_at_optimal_temperature = exchange.arrhenius_2('Vcmax', 39.0, setup_default_photosynthesis_params())
    value_at_high_temperature = exchange.arrhenius_2('Vcmax', 50.0, setup_default_photosynthesis_params())

    assert all([value_at_optimal_temperature > val for val in (value_at_low_temperature, value_at_high_temperature)])


def test_arrhenius_2_is_maximum_for_jmax_at_37_degrees_celsius():
    value_at_low_temperature = exchange.arrhenius_2('Jmax', 0.0, setup_default_photosynthesis_params())
    value_at_optimal_temperature = exchange.arrhenius_2('Jmax', 37.0, setup_default_photosynthesis_params())
    value_at_high_temperature = exchange.arrhenius_2('Jmax', 50.0, setup_default_photosynthesis_params())

    assert all([value_at_optimal_temperature > val for val in (value_at_low_temperature, value_at_high_temperature)])


def test_dhd_sensibility_decreases_as_leaf_water_potential_decreases():
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


def test_dhd_sensibility_is_not_affected_by_leaf_temperature_increases_when_leaf_water_potential_is_optimal():
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


def test_dhd_sensibility_is_affected_by_leaf_temperature_when_leaf_water_potential_is_below_its_minimum_threshold():
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


def test_dhd_sensibility_is_affected_by_leaf_temperature_only_within_restricted_bounds():
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


def test_compute_an_2par_increases_only_electron_transport_as_ppfd_increases():
    vcmax = []
    j_frac = []
    tpu_triple = []
    rd = []

    for ppfd in range(0, 2000, 10):
        res = exchange.compute_an_2par(setup_default_photosynthesis_params(), ppfd, leaf_temp=25.0)
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


def test_compute_an_2par_affects_all_photosynthetic_parameters_by_temperature():
    vcmax = []
    j_frac = []
    tpu_triple = []
    rd = []

    for temperature in range(-10, 45):
        res = exchange.compute_an_2par(setup_default_photosynthesis_params(), ppfd=1800.,
                                       leaf_temp=temperature)
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
                                              ideal_gas_cst=constants.ideal_gaz_cst) for l in arange(0.001, 0.25, 0.02)]
    assert all(x >= y for x, y in zip(gb, gb[1:]))


def test_boundary_layer_conductance_increases_as_wind_speed_increases():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=w, atm_pressure=101.3, air_temp=25.,
                                              ideal_gas_cst=constants.ideal_gaz_cst) for w in arange(0.0, 5.0, 0.5)]
    assert all(x <= y for x, y in zip(gb, gb[1:]))


def test_boundary_layer_conductance_is_weakly_dependent_on_atmospheric_pressure():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=2., atm_pressure=p, air_temp=25.,
                                              ideal_gas_cst=constants.ideal_gaz_cst) for p in arange(90.3, 102.3)]
    assert [testing.assert_almost_equal(x, y, decimal=3) for x, y in zip(gb, gb[1:])]


def test_boundary_layer_conductance_increases_as_air_temperature_increases():
    gb = [exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=2.0, atm_pressure=101.3, air_temp=t,
                                              ideal_gas_cst=constants.ideal_gaz_cst) for t in range(-10, 46)]
    assert all(x <= y for x, y in zip(gb, gb[1:]))


def test_an_gs_ci_reduces_gas_exchange_rates_as_leaf_water_potential_decreases(
        leaf_local_weather=setup_leaf_local_weather()):
    an, _, _, gs = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=25.,
            leaf_water_potential=psi,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for psi in arange(0, -3, -0.1)])

    assert all(x >= y for x, y in zip(an, an[1:]))
    assert all(x >= y for x, y in zip(gs, gs[1:]))


def test_an_gs_ci_changes_gas_exchange_rates_as_leaf_temperature_changes(leaf_local_weather=setup_leaf_local_weather()):
    an, _, _, gs = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=t,
            leaf_water_potential=0,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for t in range(-10, 46)])

    assert all(x != y for x, y in zip(an, an[1:]))
    assert all(x != y for x, y in zip(gs, gs[1:]))


def test_an_gs_ci_yields_maximum_net_photosynthesis_at_31_degrees_celsius(
        leaf_local_weather=setup_leaf_local_weather()):
    an_max, *_ = exchange.an_gs_ci(
        air_temperature=leaf_local_weather['Tac'],
        absorbed_ppfd=leaf_local_weather['PPFD'],
        relative_humidity=leaf_local_weather['hs'],
        leaf_temperature=31,
        leaf_water_potential=0,
        photo_params=setup_default_photosynthesis_params(),
        gs_params=setup_default_stomatal_conductance_params(),
        rbt=2. / 3.,
        ca=400.)

    an, *_ = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=t,
            leaf_water_potential=0,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for t in range(-10, 46)])

    assert all(x <= an_max for x in an)


def test_an_gs_ci_yields_maximum_stomatal_conductance_at_34_degrees_celsius(
        leaf_local_weather=setup_leaf_local_weather()):
    *_, gs_max = exchange.an_gs_ci(
        air_temperature=leaf_local_weather['Tac'],
        absorbed_ppfd=leaf_local_weather['PPFD'],
        relative_humidity=leaf_local_weather['hs'],
        leaf_temperature=34,
        leaf_water_potential=0,
        photo_params=setup_default_photosynthesis_params(),
        gs_params=setup_default_stomatal_conductance_params(),
        rbt=2. / 3.,
        ca=400.)

    *_, gs = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=t,
            leaf_water_potential=0,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for t in range(-10, 46)])

    assert all(x <= gs_max for x in gs)


def test_an_gs_ci_yields_more_severe_stress_when_temperature_and_water_stresses_are_combined_than_separated(
        leaf_local_weather=setup_leaf_local_weather()):
    an_t, *_, gs_t = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=t,
            leaf_water_potential=0,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for t in range(-10, 46)])

    an_t_psi, *_, gs_t_psi = zip(
        *[exchange.an_gs_ci(
            air_temperature=leaf_local_weather['Tac'],
            absorbed_ppfd=leaf_local_weather['PPFD'],
            relative_humidity=leaf_local_weather['hs'],
            leaf_temperature=t,
            leaf_water_potential=-2,
            photo_params=setup_default_photosynthesis_params(),
            gs_params=setup_default_stomatal_conductance_params(),
            rbt=2. / 3.,
            ca=400.)
            for t in range(-10, 46)])

    assert all(x <= x_t for x, x_t in zip(an_t_psi, an_t))
    assert all(x <= x_t for x, x_t in zip(gs_t_psi, gs_t))


def test_transpiration_rate_incrases_as_vapor_pressure_deficit_increases(leaf_local_weather=setup_leaf_local_weather()):
    air_temp = 25.
    leaf_temp = 25.
    atmospheric_pressure = leaf_local_weather['Pa']
    gb = exchange.boundary_layer_conductance(leaf_length=0.1, wind_speed=leaf_local_weather['u'],
                                             atm_pressure=leaf_local_weather['Pa'], air_temp=air_temp,
                                             ideal_gas_cst=constants.ideal_gaz_cst)
    es = utilities.saturated_air_vapor_pressure(leaf_temp)
    *_, gs = exchange.an_gs_ci(
        air_temperature=leaf_local_weather['Tac'],
        absorbed_ppfd=leaf_local_weather['PPFD'],
        relative_humidity=leaf_local_weather['hs'],
        leaf_temperature=34,
        leaf_water_potential=0,
        photo_params=setup_default_photosynthesis_params(),
        gs_params=setup_default_stomatal_conductance_params(),
        rbt=2. / 3.,
        ca=400.)

    transpiration = [exchange.calc_transpiration_rate(leaf_temp, ea, gs, gb, atmospheric_pressure)
                     for ea in linspace(es, 0, 10)]
    assert all(x <= y for x, y in zip(transpiration, transpiration[1:]))
