# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Gas-exchange module of HydroShoot.

This module computes net photosynthesis and stomatal conductance rates.

"""
from copy import deepcopy
from math import exp, acos, sqrt, cos, log

from hydroshoot import utilities as utils
from hydroshoot.constants import oxygen_partial_pressure as o, ideal_gaz_cst as r


# ==============================================================================
# Farquhar Parameters
# ==============================================================================

def leaf_Na(age_gdd, ppfd_10, a_n=-0.0008, b_n=3.3, a_m=6.471, b_m=56.635):
    """Computes Nitrogen content per unit leaf area.

    Args:
        age_gdd (float): [°Cd] leaf age given in cumulative degree-days temperature since budburst
        ppfd_10 (float): [umol m-2 s-1] cumulative intercepted irradiance (PPFD) over the 10 days prior to simulation
            period
        a_n (float): [gN gDM-1 °Cd-1] slope of the linear relationship between nitrogen content per unit mass
            area and leaf age
        b_n (float): [gN gDM-1] intercept of the linear relationship between nitrogen content per unit mass
            area and leaf age
        a_m (float): [gDM s umol-1] slope of the linear relationship between leaf dry matter per unit leaf area
            and absorbed ppfd over the past 10 days
        b_m (float): [gDM m-2] intercept of the linear relationship between leaf dry matter per unit leaf area
            and absorbed ppfd over the past 10 days

    Returns:
         (float): [gN m-2] nitrogen content per unit leaf area

    References:
        Prieto et al. (2012)
            A leaf gas exchange model that accounts for intra-canopy variability by considering leaf nitrogen content
            and local acclimation to radiation in grapevine (Vitis vinifera L.).
            Plant, Cell and Environment 35, 1313 - 1328

    Notes:
        Deflaut parameters' values are given for Syrah from experiments in Montpellier (France).

    """

    leaf_mass_per_area = a_m * log(max(1.e-3, ppfd_10)) + b_m
    nitrogen_content_per_leaf_mass = a_n * age_gdd + b_n
    nitrogen_content_per_leaf_area = leaf_mass_per_area * nitrogen_content_per_leaf_mass / 100.

    return max(0., nitrogen_content_per_leaf_area)


# ==============================================================================
# compute An
# ==============================================================================

def arrhenius_1(param_name, leaf_temperature, photo_params):
    """Computes the effect of temperature on the photosynthetic parameters `Tx`, `Kc`, and `Ko`.

    Args:
        param_name (str): name of the parameter to be considered, one of 'Tx', 'Kc', and 'Ko'
        leaf_temperature (float): [°C] leaf temperature
        photo_params (dict): default values at 25 °C of Farquhar's model

    Returns:
        (float): the value of the given :arg:`param_name` at the given :arg:`leaf_temperature`

    References:
        Bernacchi et al. (2003)
            In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis.
            Plant, Cell and Environment 26, 1419 –  1430

    """

    param_list = {
        'Tx': ('Tx25', 'RespT_Tx'),
        'Kc': ('Kc25', 'RespT_Kc'),
        'Ko': ('Ko25', 'RespT_Ko')}

    temp_k = utils.celsius_to_kelvin(leaf_temperature)
    param_key = param_list[param_name][1]
    shape_param, activation_energy = [photo_params[param_key][x] for x in ('c', 'deltaHa')]
    param_value = exp(shape_param - (activation_energy / (r * temp_k)))

    return param_value


def arrhenius_2(param_name, leaf_temperature, photo_params):
    """Computes the effect of temperature on the photosynthetic parameters `Vcmax`, `Jmax`, `TPUmax`, and `Rdmax`.

    Args:
        param_name (str): name of the parameter to be considered, one of `Vcmax`, `Jmax`, `TPUmax`, and `Rdmax`
        leaf_temperature (float): [°C] leaf temperature
        photo_params (dict): default values at 25 °C of Farquhar's model

    Returns:
        (float): the value of the given :arg:`param_name` at the given :arg:`leaf_temperature`

    References:
        Bernacchi et al. (2003)
            In vivo temperature response functions of parameters required to model RuBP-limited photosynthesis.
            Plant, Cell and Environment 26, 1419 –  1430

    """

    enthalpy_activation = photo_params['ds']
    enthalpy_deactivation = photo_params['dHd']

    param_list = {
        'Vcmax': ('Vcm25', 'RespT_Vcm'),
        'Jmax': ('Jm25', 'RespT_Jm'),
        'TPUmax': ('TPU25', 'RespT_TPU'),
        'Rdmax': ('Rd', 'RespT_Rd')}

    temp_k = utils.celsius_to_kelvin(leaf_temperature)
    param_key = param_list[param_name][1]
    param_value_at_25 = photo_params[param_list[param_name][0]]
    shape_param, activation_energy = [photo_params[param_key][x] for x in ('c', 'deltaHa')]
    param_value = param_value_at_25 * (exp(shape_param - (activation_energy / (r * temp_k)))) / (
            1. + exp((enthalpy_activation * temp_k - enthalpy_deactivation) / (r * temp_k)))

    return param_value


def dHd_sensibility(psi, temp, dhd_max=200.,
                    dhd_inhib_beg=195., dHd_inhib_max=190.,
                    psi_inhib_beg=-.75, psi_inhib_max=-2.,
                    temp_inhib_beg=35, temp_inhib_max=40):
    """Calculates the combined effect of irradiance and heat stress on the enthalpy of deactivation parameter.

    Args:
        psi (float): [MPa] leaf water potential
        temp (float): [°C] leaf temperature
        dhd_max (float): [KJ mol-1] maximum value of enthalpy of deactivation
        dhd_inhib_beg (float): [KJ mol-1] value of enthalpy of deactivation at the begining of photoinhibition
        dHd_inhib_max (float): [KJ mol-1] value of enthalpy of deactivation under maximum photoinhibition
        psi_inhib_beg (float): [MPa] leaf water potential at which photoinhibition begins
        psi_inhib_max (float): [MPa] leaf water potential at which photoinhibition is maximum
        temp_inhib_beg (float): [°C] leaf temperature at which photoinhibition begins
        temp_inhib_max (float): [°C] leaf temperature at which photoinhibition is maximum

    Returns:
        (float): [KJ mol-1] value of enthalpy of deactivation after considering photoinhibition

    """

    dhd_temp_effect = dhd_inhib_beg - (dhd_inhib_beg - dHd_inhib_max) * min(1.,
                                                                            max(0., (temp - temp_inhib_beg)) / float(
                                                                                temp_inhib_max - temp_inhib_beg))
    dhd_psi_effect = dhd_max - max(0., (dhd_max - dhd_temp_effect) * min(1., (psi - psi_inhib_beg) / float(
        psi_inhib_max - psi_inhib_beg)))

    return dhd_psi_effect


def compute_an_2par(params_photo, ppfd, leaf_temp):
    """Calculates the photosynthetic variables required for the analytical solution of net assimilation rate -
    stomatal conductance.

    Args:
        params_photo (dict):
        ppfd (float): [umol m-2 s-1] absorbed photosynthetic photon flux density
        leaf_temp (float): [°C] leaf temperature

    Returns:
        Vcmax (float): [umol m-2 s-1] maximum RuBP-saturated rate of carboxylation
        Kc*(1+O/Ko) (float): where
            Kc [umol mol-1] is Michaelis-Menten constant for the carboxylase
            Ko [mmol mol-1] is Michaelis-Menten constant for the oxygenase
            O [mmol mol-1] is oxygen partial pressure
        J/4. (float): where J [umol m-2 s-1] is electron transport
        2.*T (float): where T [umol mol-1] is CO2 compensation point in the absence of mitochondrial respiration
        3*TPU (float): where TPU [umol m-2 s-1] is the rate of triose phosphate transport
        -T (float): [umol mol-1] CO2 compensation point in the absence of mitochondrial respiration
        Rd (float): [umol m-2 s-1] rate of mitochondrial respiration

    References:
        Evers et al. (2010)
            Simulation of wheat growth and development based on organ-level photosynthesis and assimilate allocation.
            Journal of Experimental Botany 61, 2203 - 2216.

    """

    gamma = arrhenius_1('Tx', leaf_temp, params_photo)
    k_c = arrhenius_1('Kc', leaf_temp, params_photo)
    k_o = arrhenius_1('Ko', leaf_temp, params_photo)

    v_cmax = arrhenius_2('Vcmax', leaf_temp, params_photo)
    j_max = arrhenius_2('Jmax', leaf_temp, params_photo)
    tpu = arrhenius_2('TPUmax', leaf_temp, params_photo)
    r_d = arrhenius_2('Rdmax', leaf_temp, params_photo)

    alpha = params_photo['alpha']

    j = (alpha * ppfd) / ((1 + ((alpha ** 2 * ppfd ** 2) / (j_max ** 2))) ** 0.5)

    x1c = v_cmax
    x2c = k_c * (1 + o / k_o)

    x1j = j / 4.
    x2j = 2. * gamma

    x1t = 3 * tpu
    x2t = -gamma

    return x1c, x2c, x1j, x2j, x1t, x2t, r_d


# ==============================================================================
# Analytical solution of the An - gs - ci processes
# ==============================================================================
def fvpd_3(model, vpd, psi, psi_crit=-0.37, m0=5.278, steepness_tuzet=1.85, d0_leuning=30.):
    """Calculates the effect of water deficit on stomatal conductance.

    Args:
        model (str): stomatal conductance reduction model, one of 'misson','tuzet', 'linear' or 'vpd'
        vpd (float): [kPa] vapor pressure deficit
        psi (float): [MPa] leaf water potential
        psi_crit (float): [MPa] critical leaf water potential (the meaning varies according to the selected
            :arg:`model`)
        m0 (float): [mmol(H2O) umol-1(CO2)] slope between the stomatal conductance and assimilated CO2 rates
        steepness_tuzet (float): [MPa-1] steepness of the sigmoidal reduction function of tuzet's model
        d0_leuning (float): [kPa-1] shape factor shaping VPD's effect on the reduction function of leuning's model

    Returns:
        (float): [mmol(H2O) umol-1(CO2)] the slope between the stomatal conductance and assimilated CO2 rates after
            accounting for the effect of water deficit

    References:
        Leuning (1995)
            A critical appraisal of a combined stomatal photosynthesis model for C3 plants.
            Plant, Cell and Environment 18, 339355.
        Tuzet et al. (2003)
            A coupled model of stomatal conductance, photosynthesis and transpiration.
            Plant, Cell and Environment 26, 10971116.
        Misson et al. (2004)
            A comparison of three approaches to modeling leaf gas exchange in annually drought-stressed ponderosa pine
                forest.
            Tree Physiology 24, 529 541.

    """
    if model == 'misson':
        reduction_factor = 1. / (1. + (psi / psi_crit) ** steepness_tuzet)
    elif model == 'tuzet':
        reduction_factor = (1. + exp(steepness_tuzet * psi_crit)) / (
                1. + exp(float(steepness_tuzet * (psi_crit - psi))))
    elif model == 'linear':
        reduction_factor = 1. - min(1., float(psi / psi_crit))
    elif model == 'vpd':
        reduction_factor = 1. / (1. + vpd / float(d0_leuning))
    else:
        raise ValueError("The 'model' argument must be one of the following ('misson','tuzet', 'linear' or 'vpd').")
    return m0 * reduction_factor


def mesophyll_conductance(leaf_temperature, gm25=0.1025, activation_energy=49600., deactivation_energy=437400.,
                          entropy=1400.):
    """Calculates mesophyll conductance to CO2.

    Args:
        leaf_temperature (float): [°C] leaf temperature
        gm25 (float): [mol m-2 s-1 bar-1] mesophyll conductance for CO2 at 25 °C
        activation_energy (float): [J mol-1] activation energy
        deactivation_energy (float): [J mol-1] deactivation energy
        entropy (float): [J K-1 mol -1]

    Returns:
        (float): [mol m-2 s-1] mesophyll conductance to CO2

    """

    R2 = r * 1.e3  # Ideal gaz constant in [J K-1 mol-1]

    return gm25 * exp((1. / 298. - 1. / (leaf_temperature + 273.)) * activation_energy / R2) * (
            1. + exp(entropy / R2 - deactivation_energy / 298. / R2)) / (
                   1. + exp(entropy / R2 - 1. / (leaf_temperature + 273.) * deactivation_energy / R2))


def boundary_layer_conductance(leaf_length, wind_speed, atm_pressure, air_temp, ideal_gas_cst):
    """Computes boundary layer conductance to water vapor.

    Args:
        leaf_length (float): [m] leaf length
        wind_speed (float): [m s-1] local wind speed
        atm_pressure (float): [kPa] atmospheric pressure
        air_temp (float): [°C] air temperature
        ideal_gas_cst (float): [kJ K-1 mol-1] Ideal gaz constant

    Returns:
        (float): [mol m-2 s-1] boundary layer conductance

    References:
        Nobel P. 2005.
            Temperature and energy budgets.
            In Nobel S, eds. Physicochemical and Environmental Plant Physiology.
            Elsevier Academic Press, 307–350.

    """

    l_w = leaf_length / 100. * 0.72  # effective length in the downwind direction [m]
    d_bl = 4. * (l_w / max(1.e-3, wind_speed)) ** 0.5 / 1000.  # Boundary layer thickness in [m] (Nobel, 2009 pp.337)
    dj0 = 2.13 * 1.e-5  # [m2 s-1] at P=1. atm and t=0. °C (Nobel, pp.545)
    dj = dj0 * (101.3 / atm_pressure) * ((air_temp + 273.15) / 273.15) ** 1.8  # (Nobel, eq.8.8, pp.379)
    gb = dj * (atm_pressure * 1.e-3) / (
            (ideal_gas_cst * 1.e-3) * (air_temp + 273.15) * d_bl)  # [mol m-2 s-1] (Nobel, 2009 pp.337)

    return gb


def compute_amono_analytic(x1, x2, leaf_temperature, vpd, gammax, rd, psi, model='misson',
                           g0=0.019, rbt=2. / 3., ca=400., m0=5.278, psi0=-0.1,
                           d0_leuning=30., steepness_tuzet=1.85):
    """Computes the gross CO2 assimilation rate analytically.

    Args:
        x1 (float): Vcmax or J/4.
        x2 (float): KmC(1+O/KmO) or 2 gamma_asterisk
        leaf_temperature (float): [°C] leaf temperature
        vpd (float): [kPa] leaf-to-air vapor pressure deficit
        gammax (float): [ubar] CO2 compensation point
        rd (float): [umol m-2 s-1] mitochondrial respiration rate in the light
        psi (float): [MPa] bulk water potential of the leaf
        model (str): stomatal conductance reduction model, one of 'misson','tuzet', 'linear' or 'vpd'
        g0 (float): [umol m-2 s-1] residual stomatal conductance to CO2
        rbt (float): [m2 s ubar umol-1] combined turbulance and boundary layer resistance to CO2
        ca (float): [ubar] CO2 partial pressure of the air
        m0 (float): [umol mmol-1] maximum slope An/gs (absence of water deficit), see :func:`fvpd_3`
        psi0 (float): [MPa] critical thershold for water potential, see :func:`fvpd_3`
        d0_leuning (float): [kPa-1] shape parameter, see :func:`fvpd_3`
        steepness_tuzet (float): [MPa-1] shape parameter, see :func:`fvpd_3`

    Returns:
        (float): [umolCO2 m-2 s-1] gross CO2 assimilation rate

    References:
        Evers et al. 2010.
            Simulation of wheat growth and development based on organ-level photosynthesis and assimilate allocation.
            Journal of Experimental Botany 61, 2203 – 2216.

    Notes:
        A is called 'mono' since it is calculated for either of RuBisco- or e transport-limited phases.
        The net A (An) is calculated using the function 'compute_an_analytic'

    """

    f_vpd = fvpd_3(model, vpd, psi, psi_crit=psi0, m0=m0, steepness_tuzet=steepness_tuzet, d0_leuning=d0_leuning)

    g_m = mesophyll_conductance(leaf_temperature)
    cube_a = g0 * (x2 + gammax) + (g0 / g_m + f_vpd) * (x1 - rd)
    cube_b = utils.cmol2cpa(leaf_temperature, ca) * (x1 - rd) - gammax * x1 - rd * x2
    cube_c = utils.cmol2cpa(leaf_temperature, ca) + x2 + (1. / g_m + rbt) * (x1 - rd)
    cube_d = x2 + gammax + (x1 - rd) / g_m
    cube_e = 1. / g_m + (g0 / g_m + f_vpd) * (1. / g_m + rbt)
    cube_p = -(cube_d + (x1 - rd) / g_m + cube_a * (1. / g_m + rbt) + (g0 / g_m + f_vpd) * cube_c) / cube_e
    cube_q = (cube_d * (x1 - rd) + cube_a * cube_c + (g0 / g_m + f_vpd) * cube_b) / cube_e
    cube_r = -(cube_a * cube_b / cube_e)
    cube_Q = max(1.e-32, (cube_p ** 2. - 3. * cube_q) / 9.)
    cube_R = (2. * cube_p ** 3. - 9. * cube_p * cube_q + 27. * cube_r) / 54.
    cube_xi = acos(max(-1., min(1., cube_R / sqrt(cube_Q ** 3.))))

    a_mono = -2. * sqrt(cube_Q) * cos(cube_xi / 3.) - cube_p / 3.

    return a_mono


def compute_an_analytic(leaf_temperature, vpd, x1c, x2c, x1j, x2j, x1t, x2t, rd, psi,
                        model='misson', g0=0.019, rbt=2. / 3., ca=400., m0=5.278,
                        psi0=-0.1, d0_leuning=30., steepness_tuzet=1.85):
    """Computes the gross CO2 assimilation rate analytically.

    Args:
        leaf_temperature (float): [°C] leaf temperature
        vpd (float): [kPa] leaf-to-air vapor pressure deficit
        x1c (float): [umol m-2 s-1] Rubisco-limiting CO2 assimilation parameter Vcmax
        x2c (float): [umol m-2 s-1] Rubisco-limiting CO2 assimilation parameter Kc*(1+O/Ko)
        x1j (float): [umol m-2 s-1] electron transport-limiting CO2 assimilation parameter J/4
        x2j (float): [umol m-2 s-1] electron transport-limiting CO2 assimilation parameter 2*gammax
        x1t (float): [umol m-2 s-1] triose phosphate export rates parameter 3*TPU
        x2t (float): [umol m-2 s-1] triose phosphate export rates parameter -gamma
        rd (float): [umol m-2 s-1] mitochondrial respiration rate in the light
        psi (float): [MPa] bulk water potential of the leaf
        model (str): stomatal conductance reduction model, one of 'misson','tuzet', 'linear' or 'vpd'
        g0 (float): [umol m-2 s-1] residual stomatal conductance to CO2
        rbt (float): [m2 s ubar umol-1] combined turbulance and boundary layer resistance to CO2
        ca (float): [ubar] CO2 partial pressure of the air
        m0 (float): [umol mmol-1] maximum slope An/gs (absence of water deficit), see :func:`fvpd_3`
        psi0 (float): [MPa] critical thershold for water potential, see :func:`fvpd_3`
        d0_leuning (float): [kPa-1] shape parameter, see :func:`fvpd_3`
        steepness_tuzet (float): [MPa-1] shape parameter, see :func:`fvpd_3`

    Returns:
        (float): [umolCO2 m-2 s-1] net CO2 assimilation rate
        (float): [ubar] CO2 chloroplast partial pressure
        (float): [ubar] intercullar CO2 partial pressure
        (float): [mol m-2 s-1] stomatal conductance to water vapor

    """

    f_vpd = fvpd_3(model, vpd, psi, psi_crit=psi0, m0=m0, steepness_tuzet=steepness_tuzet, d0_leuning=d0_leuning)

    gammax = x2j / 2.

    ci_asterisk = gammax - rd / mesophyll_conductance(leaf_temperature)

    a_c = compute_amono_analytic(x1c, x2c, leaf_temperature, vpd, gammax, rd, psi, model, g0, rbt, ca, m0, psi0,
                                 d0_leuning, steepness_tuzet)
    a_j = compute_amono_analytic(x1j, x2j, leaf_temperature, vpd, gammax, rd, psi, model, g0, rbt, ca, m0, psi0,
                                 d0_leuning, steepness_tuzet)
    a_t = compute_amono_analytic(x1t, x2t, leaf_temperature, vpd, gammax, rd, psi, model, g0, rbt, ca, m0, psi0,
                                 d0_leuning, steepness_tuzet)

    a_n = min(a_c, a_j, a_t)

    # chlorophyll partial pressure [ubar]
    if min(a_c, a_j, a_t) == a_c:
        c_c = (gammax * x1c + (a_c + rd) * x2c) / (x1c - a_c - rd)
    elif min(a_c, a_j, a_t) == a_j:
        c_c = (gammax * x1j + (a_j + rd) * x2j) / (x1j - a_j - rd)
    else:
        c_c = (gammax * x1t + (a_t + rd) * x2t) / (x1t - a_t - rd)

    # inter-cellular partial pressure [ubar]
    c_i = c_c + a_n / mesophyll_conductance(leaf_temperature)

    # stomatal conductance to CO2 [umol m-2 s-1 ubar-1]
    gs = g0 / 1.6 + (a_n + rd) * f_vpd / (c_i - ci_asterisk)

    # stomatal conductance to water vapor [umol m-2 s-1 ubar-1]
    gsw = gs * 1.6

    return a_n, utils.cpa2cmol(leaf_temperature, c_c), utils.cpa2cmol(leaf_temperature, c_i), gsw


def an_gs_ci(absorbed_ppfd, vapor_pressure_deficit, leaf_temperature, leaf_water_potential,
             photo_params, gs_params, rbt=2. / 3., ca=400.):
    """Computes simultaneously the net CO2 assimilation rate (An), stomatal conductance to water vapor (gs), and
    inter-cellular CO2 concentration (Ci), by a_n analytic scheme.

    Args:
        absorbed_ppfd (float): [umol m-2 s-1] absorbed Photosynthetic Photon Flux Density
        vapor_pressure_deficit (float): (kPa) leaf-to-air vapor pressure deficit
        leaf_temperature (float): [°C] leaf temperature
        leaf_water_potential (float): [MPa] bulk water potential of the leaf
        photo_params (dict): default values at 25 °C of Farquhar's model
        gs_params (dict): stomatal conductance params (keys must be: 'model', 'g0', 'm0', 'psi0', 'D0', 'n')
        rbt (float): [m2 s ubar umol-1] combined turbulance and boundary layer resistance to CO2
        ca (float): [ubar] CO2 partial pressure of the air

    Returns: (tuple)
        (float): [umolCO2 m-2 s-1] net CO2 assimilation rate
        (float): [ubar] CO2 chloroplast partial pressure
        (float): [ubar] intercullar CO2 partial pressure
        (float): [mol m-2 s-1] stomatal conductance to water vapor

    """
    x1c, x2c, x1j, x2j, x1t, x2t, r_d = compute_an_2par(
        params_photo=photo_params,
        ppfd=max(1.e-6, absorbed_ppfd),  # To avoid numerical instability
        leaf_temp=leaf_temperature)

    a_n, c_c, c_i, gs = compute_an_analytic(
        leaf_temperature=leaf_temperature,
        vpd=vapor_pressure_deficit,
        psi=leaf_water_potential,
        x1c=x1c,
        x2c=x2c,
        x1j=x1j,
        x2j=x2j,
        x1t=x1t,
        x2t=x2t,
        rd=r_d,
        model=gs_params['model'],
        g0=gs_params['g0'],
        m0=gs_params['m0'],
        psi0=gs_params['psi0'],
        d0_leuning=gs_params['D0'],
        steepness_tuzet=gs_params['n'],
        rbt=rbt,
        ca=ca)

    return a_n, c_c, c_i, gs


def calc_transpiration_rate(vpd, gs, gb, atm_pressure=101.3):
    """Computes transpiration rate per unit leaf surface area.

    Args:
        vpd (float): [kPa] leaf-to-air vapor pressure deficit
        gs (float): [mol m-2 s-1] stomatal conductance to water vapor
        gb (float): [mol m-2 s-1] boundary layer conductance to water vapor
        atm_pressure (float): [kPa] atmospheric pressure

    Returns:
        (float): [mol m-2leaf s-1] transpiration rate

    Warnings:
        gb is calculated for both sides of a flat leaf.

    """

    gv = 1. / ((2. / gb) + (1. / gs))
    return gv * (vpd / atm_pressure)


def calc_gas_exchange_rates(leaf_water_potential, leaf_temperature, carboxylation_max, electron_transport,
                            triose_photsphate_transport, mitochondreal_respiration,
                            entalpy_deactivation, absorbed_ppfd, leaf_length, wind_speed, photo_params, gs_params,
                            air_temperature, vapor_pressure_deficit, air_co2, atmospheric_pressure, rbt=2. / 3.,
                            leaf_ids=None):
    """Computes gas exchange fluxes at the leaf scale analytically.

    Args:
        leaf_water_potential (dict): [MPa] leaf water potential
        leaf_temperature (dict): [°C] leaf temperature
        carboxylation_max (dict): [umol m-2 s-1] maximum RuBP-saturated rate of carboxylation
        electron_transport (dict): [umol m-2 s-1] electron transport
        triose_photsphate_transport (dict): [umol m-2 s-1] triose phosphate transport rate
        mitochondreal_respiration (dict): [umol m-2 s-1] mitochondrial respiration rate
        entalpy_deactivation (dict): [KJ mol-1] enthalpy of deactivation
        absorbed_ppfd (dict): [umol m-2 s-1] absorbed photosynthetic photon flux density
        leaf_length (dict): [m] leaf length
        wind_speed (dict): [m s-1] wind speed at the leaf level
        leaf_ids (list): leaf ids (default None)
        photo_params (dict): values at 25 °C of Farquhar's model
        gs_params (dict): parameters of the stomatal conductance model (model, g0, m0, psi0, D0, n)
        air_temperature (dict): [°C] air temperature
        vapor_pressure_deficit (float): [kPa] leaf-to-air vapor pressure deficit
        air_co2 (float): [ppm] air CO2 concentration
        atmospheric_pressure (float): [kPa] atmospheric pressure
        rbt (float): [m2 s ubar umol-1] the combined turbulance and boundary layer resistance to CO2 transport

    References:
        Evers et al. 2010.
            Simulation of wheat growth and development based on organ-level photosynthesis and assimilate allocation.
            Journal of Experimental Botany 61, 2203 – 2216.

    Returns:
        (dict): [umol m-2 s-1] the net CO2 assimilation
        (dict): [umol mol] intercellular CO2 concentration
        (dict): [mol m-2 s-1] stomatal conductance to water vapor
        (dict): [mol m-2 s-1] boundary layer conductance to water vapor
        (dict): [mol m-2leaf s-1] transpiration per unit leaf surface area

    """
    net_photosynthesis_rate = {}
    intercellular_carbon_content = {}
    stomatal_conductance_rate = {}
    boundary_layer_conductance_rate = {}
    transpiration_rate = {}

    for vid in (leaf_ids if leaf_ids is not None else leaf_length.keys()):
        water_potential_leaf = leaf_water_potential[vid]
        temperature_leaf = leaf_temperature[vid]
        vpd_leaf = vapor_pressure_deficit[vid]
        entalpy_deactivation[vid] = dHd_sensibility(
            psi=water_potential_leaf,
            temp=temperature_leaf,
            dhd_max=photo_params['dHd'],
            dhd_inhib_beg=photo_params['photo_inhibition']['dhd_inhib_beg'],
            dHd_inhib_max=photo_params['photo_inhibition']['dHd_inhib_max'],
            psi_inhib_beg=photo_params['photo_inhibition']['psi_inhib_beg'],
            psi_inhib_max=photo_params['photo_inhibition']['psi_inhib_max'],
            temp_inhib_beg=photo_params['photo_inhibition']['temp_inhib_beg'],
            temp_inhib_max=photo_params['photo_inhibition']['temp_inhib_max'])

        leaf_par_photo = deepcopy(photo_params)
        leaf_par_photo.update(dict(
            Vcm25=carboxylation_max[vid],
            Jm25=electron_transport[vid],
            TPU25=triose_photsphate_transport[vid],
            Rd=mitochondreal_respiration[vid],
            dHd=entalpy_deactivation[vid]))

        a_n, c_c, c_i, gs = an_gs_ci(
            absorbed_ppfd=absorbed_ppfd[vid],
            vapor_pressure_deficit=vpd_leaf,
            leaf_temperature=temperature_leaf,
            photo_params=leaf_par_photo,
            gs_params=gs_params,
            leaf_water_potential=water_potential_leaf,
            rbt=rbt,
            ca=air_co2)

        gb = boundary_layer_conductance(
            leaf_length=leaf_length[vid],
            wind_speed=wind_speed[vid],
            atm_pressure=atmospheric_pressure,
            air_temp=air_temperature[vid],
            ideal_gas_cst=r)

        # Transpiration
        e = calc_transpiration_rate(
            vpd=vpd_leaf,
            gs=gs,
            gb=gb,
            atm_pressure=atmospheric_pressure)

        net_photosynthesis_rate.update({vid: a_n})
        intercellular_carbon_content.update({vid: c_i})
        stomatal_conductance_rate.update({vid: gs})
        boundary_layer_conductance_rate.update({vid: gb})
        transpiration_rate.update({vid: max(0., e)})

    return (net_photosynthesis_rate, intercellular_carbon_content, stomatal_conductance_rate,
            boundary_layer_conductance_rate, transpiration_rate)


def set_gas_exchange_rates(g, photo_params, gs_params, relative_humidity, air_co2,
                           atmospheric_pressure, E_type2, leaf_ids, rbt=2. / 3.):
    """Sets gas exchange fluxes at the leaf scale analytically.

    Args:
        g: a multiscale tree graph object
        photo_params (dict): values at 25 °C of Farquhar's model
        gs_params (dict): parameters of the stomatal conductance model (model, g0, m0, psi0, D0, n)
        relative_humidity (float): (%) air relative humidity (between 0 and 100 for dry and saturated air, respectively)
        air_co2 (float): [ppm] air CO2 concentration
        atmospheric_pressure (float): [kPa] atmospheric pressure
        E_type2 (str): one of 'Ei' (intercepted irradiance) or 'Eabs' (absorbed irradiance)
        leaf_ids (list of int): leaf ids
        rbt (float): [m2 s ubar umol-1] the combined turbulance and boundary layer resistance to CO2 transport

    References:
        Evers et al. 2010.
            Simulation of wheat growth and development based on organ-level photosynthesis and assimilate allocation.
            Journal of Experimental Botany 61, 2203 – 2216.

    Notes:
        This function adds to each leaf the following properties:
            An (float): [umol m-2 s-1] the net CO2 assimilation
            Ci (float): [umol mol] intercellular CO2 concentration
            gs (float): [mol m-2 s-1] stomatal conductance to water vapor
            gb (float): [mol m-2 s-1] boundary layer conductance to water vapor
            E (float): [mol m-2leaf s-1] transpiration per unit leaf surface area

    """
    all_leaf_water_potential = g.property('psi_head')
    all_leaf_temperature = g.property('Tlc')
    all_vcm25 = g.property('Vcm25')
    all_jm25 = g.property('Jm25')
    all_tpu25 = g.property('TPU25')
    all_rd = g.property('Rd')
    all_dhd = g.property('dHd')
    all_absorbed_ppfd = g.property(E_type2)
    all_leaf_length = g.property('Length')
    all_wind_speed = g.property('u')
    all_air_temperature = g.property('Tac')
    all_vapor_pressure_deficit = g.property('vpd')

    (g.properties()['An'], g.properties()['Ci'], g.properties()['gs'], g.properties()['gb'],
     g.properties()['E']) = calc_gas_exchange_rates(
        leaf_water_potential=all_leaf_water_potential,
        leaf_temperature=all_leaf_temperature,
        carboxylation_max=all_vcm25,
        electron_transport=all_jm25,
        triose_photsphate_transport=all_tpu25,
        mitochondreal_respiration=all_rd,
        entalpy_deactivation=all_dhd,
        absorbed_ppfd=all_absorbed_ppfd,
        leaf_length=all_leaf_length,
        wind_speed=all_wind_speed,
        leaf_ids=leaf_ids,
        photo_params=photo_params,
        gs_params=gs_params,
        air_temperature=all_air_temperature,
        vapor_pressure_deficit=all_vapor_pressure_deficit,
        air_co2=air_co2,
        atmospheric_pressure=atmospheric_pressure,
        rbt=rbt)
    pass
