# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Gas-exchange module of HydroShoot.

This module computes net photosynthesis and stomatal conductance rates.

"""
from copy import deepcopy
from scipy import exp, arccos, sqrt, cos, log

from hydroshoot import utilities as utils

# Constants
O = 210 # Oxygen partial pressure [mmol mol-1]
R = 0.0083144598 # Ideal gaz constant [kJ K-1 mol-1]


#==============================================================================
# Farquhar Parameters
#==============================================================================

# Initial values (before accountinf for Na)
def par_photo_default(Vcm25= 89.0, Jm25 = 143.0, cRd = 0.008, TPU25 = 10.0,
       Kc25 = 404.9, Ko25 = 278.4, Tx25 = 42.75,
       alpha =  [0.2, 0.2, 0.19, 0.19, 0.14, 0.12], alpha_T_limit = [15, 20, 25, 30, 34, 50],
       a1 = 0.98, a2 = 0.98, a3 = 0.98, ds = 0.635, dHd = 200,
       RespT_Kc = {'model':'Ahrenius', 'c':38.05, 'deltaHa':79.43},
       RespT_Ko = {'model':'Ahrenius', 'c':20.30, 'deltaHa':36.38},
       RespT_Vcm = {'model':'Ahrenius', 'c':26.35, 'deltaHa':65.33},
       RespT_Jm = {'model':'Ahrenius', 'c':17.57, 'deltaHa':43.54},
       RespT_TPU = {'model':'Ahrenius', 'c':21.46, 'deltaHa':53.1},
       RespT_Rd = {'model':'Ahrenius', 'c':18.72, 'deltaHa':46.39},
       RespT_Tx = {'model':'Ahrenius', 'c':19.02, 'deltaHa':37.83}):
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
   - **alpha**, **alpha_T_limit**, **a1**, **a2**, **a3**: parameters regulating electron transport for sunlit and shaded leaves (**not used**)
   - **c**: empirical parameter defining the temperature response curves of each of Kc, Ko, Vcm, Jm, TPU, Rd and Tx
   - **deltaHa**: Activation energy of the Arrhenius functions [kJ molCO2-1] 
   - **ds**: float, enthalpy of activation [KJ mol-1]
   - **dHd**: float, enthalpy of deactivation [KJ mol-1]
   """

    par_photodef = {}
    par_photodef['Vcm25']= Vcm25
    par_photodef['Jm25'] = Jm25
    par_photodef['Rd'] = cRd*Vcm25
    par_photodef['TPU25'] = TPU25
    par_photodef['Kc25'] = Kc25
    par_photodef['Ko25'] = Ko25
    par_photodef['Tx25'] = Tx25
    par_photodef['alpha'] = alpha
    par_photodef['alpha_T_limit'] = alpha_T_limit
    par_photodef['a1'] = a1
    par_photodef['a2'] = a2
    par_photodef['a3'] = a3
    par_photodef['ds'] = ds
    par_photodef['dHd'] = dHd
    par_photodef['RespT_Kc'] = RespT_Kc
    par_photodef['RespT_Ko'] = RespT_Ko
    par_photodef['RespT_Vcm'] = RespT_Vcm
    par_photodef['RespT_Jm'] = RespT_Jm
    par_photodef['RespT_TPU'] = RespT_TPU
    par_photodef['RespT_Rd'] = RespT_Rd
    par_photodef['RespT_Tx'] = RespT_Tx

    return par_photodef


def par_25_N_dict(Vcm25_N = (34.02, -3.13), Jm25_N = (78.27, -17.3),
             Rd_N = (0.42, -0.01), TPU25_N = (6.24, -1.92)):
    """
    Generates a dictionary containing the (slope, intercept) values of
    the linear relationship between photosynthetic capacity parameters
    (Vcmax, Jmax, TPU, Rd) and surface-based leaf Nitrogen content,
    according to Prieto et al. (2012, doi: 10.1111/j.1365-3040.2012.02491.x)

    :Parameters:
    - **The first value of each tuple**: float, the **slope** (parameter/Na) [umol g-1 s-1]
    - **The second value of each tuple**: float, the **intercept** [umol m-2 s-1]
    """
    par_dict = {
    'Vcm25_N' : Vcm25_N,
    'Jm25_N' : Jm25_N,
    'Rd_N' : Rd_N,
    'TPU25_N' : TPU25_N}

    return par_dict


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


#==============================================================================
# compute An
#==============================================================================

def arrhenius_1(param_name, leaf_temperature, photo_params):
    """Computes the effect of temperature on the photosynthetic parameters `Tx`, `Kc`, and `Ko`.

    Args:
        param_name (str): name of the parameter to be considered, one of 'Tx', 'Kc', and 'Ko'
        leaf_temperature (float): [°C] leaf temperature
        photo_params (dict): default values at 25 °C of Farquhar's model (cf. :func:`par_photo_default`)

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
    param_value = exp(shape_param - (activation_energy / (R * temp_k)))

    return param_value


def arrhenius_2(param_name, leaf_temperature, photo_params):
    """Computes the effect of temperature on the photosynthetic parameters `Vcmax`, `Jmax`, `TPUmax`, and `Rdmax`.

    Args:
        param_name (str): name of the parameter to be considered, one of `Vcmax`, `Jmax`, `TPUmax`, and `Rdmax`
        leaf_temperature (float): [°C] leaf temperature
        photo_params (dict): default values at 25 °C of Farquhar's model (cf. :func:`par_photo_default`)

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
    param_value = param_value_at_25 * (exp(shape_param - (activation_energy / (R * temp_k)))) / (
            1. + exp((enthalpy_activation * temp_k - enthalpy_deactivation) / (R * temp_k)))

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
        -T (float): [umol m-2 s-1] rate of triose phosphate transport
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

    alpha = .24

    j = (alpha * ppfd) / ((1 + ((alpha ** 2 * ppfd ** 2) / (j_max ** 2))) ** 0.5)

    x1c = v_cmax
    x2c = k_c * (1 + O / k_o)

    x1j = j / 4.
    x2j = 2. * gamma

    x1t = 3 * tpu
    x2t = -gamma

    return x1c, x2c, x1j, x2j, x1t, x2t, r_d


#==============================================================================
# Analytical solution of the An - gs - ci processes
#==============================================================================
def fvpd_3(model, vpd, psi, psi_crit=-0.37, c1=5.278, c2=1.85, D0=30.):
    """
    """
    if model =='misson':
        m = 1. / (1. + (psi/psi_crit)** c2)
    elif model == 'tuzet':
        m = (1. + exp(c2 * psi_crit)) / (1. + exp(float(c2 * (psi_crit - psi))))
    elif model == 'linear':
        m = 1 - min(1,float(psi/psi_crit))
    elif model == 'vpd':
        m = 1. / (1. + vpd / float(D0))
    else:
        raise ValueError ("The 'model' argument must be one of the following ('misson','tuzet', 'linear').")
    return c1 * m


def g0_sensibility(psi, psi_crit=-0.1, n=3):
    """
    """
    return 1. / (1. + (psi/psi_crit)**n)


def gm(temp, gm25=0.1025, EA_gm=49600., DEA_gm=437400., Sgm=1400.):
    """
    Returns gm, the mesophyll conductance to CO2 diffusion [umol m-2 s-1 ubar-1].

    :Parameters:
    - **temp**: float, leaf temperature [degreeC]
    - **gm25**: Mesophyll conductance for CO2 diffusion [mol m-2 leaf s-1 bar-1]. Default value from (Bernacci 2002 for tobaco)
    - **EA_gm**: the activation energy for parameter X [J mol-1]; Yin 2009; Bernacci 2002 plant physiology
    - **DEA_gm**: the deactivation energy, resource?
    - **Sgm**: an entropy term [J K-1 mol -1]
    """

    R2 = R * 1.e3 #Ideal gaz constant in [J K-1 mol-1]
    gm = gm25*exp((1./298.-1./(temp+273.))*EA_gm/R2)*(1.+exp(Sgm/R2-DEA_gm/298./R2))/(1.+exp(Sgm/R2-1./(temp+273.)*DEA_gm/R2))

    return gm


def boundary_layer_conductance(leaf_length, wind_speed, atm_pressure, air_temp, ideal_gas_cst):
    """Computes boundary layer conductance of the leaf.

    :params:
    - **leaf_length**: [m] leaf length
    - **wind_speed**: [m s-1] local wind speed
    - **atm_pressure**: [kPa] atmospheric pressure
    - **air_temp**: [°C] air temperature
    - **ideal_gas_cst**: [kJ K-1 mol-1] Ideal gaz constant

    :returns:
    - [mol m-2 s-1] boundary layer conductance

    """
    l_w = leaf_length / 100. * 0.72  # effective length in the downwind direction [m]
    d_bl = 4. * (l_w / max(1.e-3, wind_speed)) ** 0.5 / 1000.  # Boundary layer thickness in [m] (Nobel, 2009 pp.337)
    Dj0 = 2.13 * 1.e-5  # [m2 s-1] at P=1. atm and t=0. °C (Nobel, pp.545)
    Dj = Dj0 * (101.3 / atm_pressure) * ((air_temp + 273.15) / 273.15) ** 1.8  # (Nobel, eq.8.8, pp.379)
    gb = Dj * (atm_pressure * 1.e-3) / ((ideal_gas_cst * 1.e-3) * (air_temp + 273.15) * d_bl)  # [mol m-2 s-1] (Nobel, 2009 pp.337)

    return gb


def compute_amono_analytic(x1, x2, temp, vpd, gammax, Rd, psi, model='misson',
                           g0=0.019, rbt=2./3., ca=400., m0=5.278, psi0=-0.1,
                           D0=30., n=1.85):
    """
    Returns Amono, the gross CO2 assimilation rate [umolCO2 m-2 s-1] according to the analytical solution of Evers et al. (2010) jxb.
    
    A is called 'mono' since it is calculated for either of RuBisco- or e transport-limited phases. The net A (An) is calculated using the function 'compute_an_analytic'
    
    :Parameters:
    - **x1** and **x2**: float, are Farquhar's An parameters (x1= Vcmax or J/4.; x2= KmC(1+O/KmO) or 2 gamma_asterisk)
    - **temp**: float, leaf temperature [degreeC]
    - **vpd**: float, vapor pressure deficit [kPa]
    - **gammax**: float, CO2 compensation point [ubar]
    - **Rd**: float, the mitochondrial respiration rate in the light [umol m-2 s-1]
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters
    """

    f_VPD = fvpd_3(model, vpd, psi, psi_crit=psi0, c1=m0, c2=n, D0=D0)
#    if psi < -1.5: f_VPD = 0

    cube_a = g0*(x2+gammax)+(g0/gm(temp)+f_VPD)*(x1-Rd)
    cube_b = utils.cmol2cpa(temp, ca) * (x1 - Rd) - gammax * x1 - Rd * x2
    cube_c = utils.cmol2cpa(temp, ca) + x2 + (1. / gm(temp) + rbt) * (x1 - Rd)
    cube_d = x2+gammax+(x1-Rd)/gm(temp)
    cube_e = 1./gm(temp)+(g0/gm(temp)+f_VPD)*(1./gm(temp)+rbt)
    cube_p = -(cube_d+(x1-Rd)/gm(temp)+cube_a*(1./gm(temp)+rbt)+(g0/gm(temp)+f_VPD)*cube_c)/cube_e
    cube_q = (cube_d*(x1-Rd)+cube_a*cube_c+(g0/gm(temp)+f_VPD)*cube_b)/cube_e
    cube_r = -(cube_a*cube_b/cube_e)
    cube_Q = (cube_p**2.-3.*cube_q)/9.
    cube_R = (2.*cube_p**3.-9.*cube_p*cube_q+27.*cube_r)/54.
    cube_xi = arccos(max(-1.,min(1., cube_R/sqrt(cube_Q**3.))))

    Amono = -2.*sqrt(cube_Q)*cos(cube_xi/3.)-cube_p/3.

    return Amono


def compute_an_analytic(temp, vpd, x1c, x2c, x1j, x2j, x1t, x2t, Rd, psi,
                        model='misson', g0=0.019, rbt=2./3., ca=400., m0=5.278,
                        psi0=-0.1, D0=30., n=1.85):
    """
    Returns An, the net CO2 assimilation rate [umolCO2 m-2 s-1] according to the analytical solution of Evers et al. (2010) jxb.
    
    :Parameters:
    - **temp**: float, leaf temperature [degreeC]
    - **vpd**: float, vapor pressure deficit [kPa]
    - **x1c** and **x2c** are Rubisco-limiting CO2 assimilation of Farquhar model (Vcmax and Kc*(1+O/Ko), resp.)
    - **x1j** and **x2j** are electron transport-limiting CO2 assimilation of Farquhar model (J/4 and 2*gammax, resp.)
    - **x1t** and **x1t** are triose phosphate export rates of Farquhar model (3*TPU and -gamma, resp.)
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters
    
    :Return:
    - **An**: CO2 net assimilation rate [umolCO2 m-2 s-1]
    - **Cc**: CO2 chloroplast partial pressure [ubar]
    - **Ci**: intercullar CO2 partial pressure [ubar]
    - **gsw**: stomatal conductance for water [mol m-2 s-1]
    """

    f_VPD = fvpd_3(model, vpd, psi, psi_crit=psi0, c1=m0, c2=n, D0=D0)

    gammax = x2j/2.
    ci_asterisk = gammax-Rd/gm(temp)

    Ac = compute_amono_analytic(x1c,x2c,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    Aj = compute_amono_analytic(x1j,x2j,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    At = compute_amono_analytic(x1t,x2t,temp,vpd,gammax,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)
    
    An = min(Ac,Aj,At)

    # chlorophyll partial pressure [ubar]
    if min(Ac,Aj,At) == Ac:
        CC=(gammax*x1c+(Ac+Rd)*x2c)/(x1c-Ac-Rd)
    elif min(Ac,Aj,At) == Aj:
        CC=(gammax*x1j+(Aj+Rd)*x2j)/(x1j-Aj-Rd)
    else:
        CC=(gammax*x1t+(At+Rd)*x2t)/(x1t-At-Rd)

    # Intercellular partial pressure [ubar]
    CI = CC + An/gm(temp);  #intercelluar co2 concentration

    # stomatal conductance for CO2 [umol m-2 s-1 ubar-1]
    GS =g0/1.6+(An + Rd)*f_VPD/(CI-ci_asterisk)
    # stomatal conductance for water [umol m-2 s-1 ubar-1]
    GSW = GS*1.6

    return An, utils.cpa2cmol(temp, CC), utils.cpa2cmol(temp, CI), GSW


def an_gs_ci(par_photo, meteo_leaf, psi, Tlc, model='misson', g0=0.019,rbt=2./3.,
           ca=400., m0=5.278, psi0=-0.1, D0=30., n=1.85):
    """
    Estimates simultaneously the values of net CO2 assimilation rate (An) and
    intercellular CO2 concentration (Ci), by an analytic scheme.

    :Parameters:
    - **par_photo**: a dictionary of Farquhar's model parameters
    - **meteo_leaf**: list of local (leaf-scale) meteorological data
    - **psi**: float, bulk water potential of the leaf [MPa]
    - **Tlc**: float, leaf temperature [degreeC]
    - **model**: string, one of 'misson', 'tuzet', 'linear', the response function of stomata to water status
    - **g0**: float, the residual stomatal conductance for CO2 at the light compensation point [umol m-2 s-1]
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    - **m0**: float, the maximum slope An/gs (absence of water deficit) [umol mmol-1]
    - **psi0**, : float, a critical thershold for water potential [MPa]
    - **D0** and **n**: float, shape parameters

    :Returns:
    - **An**: net CO2 assimilation rate [umol m-2 s-1]
    - **Cc**: chloroplast CO2 concentration [umol mol-1]
    - **Cc**: intercellular CO2 concentration [umol mol-1]
    - **gs**: stomatal conductances to water [mol m-2 s-1]
    """

    Tac = meteo_leaf['Tac']
    PPFD = meteo_leaf['PPFD'] 
    hs = meteo_leaf['hs']
    
    PPFD = max(1.e-6, PPFD) # To avoid numerical instability at PPFD=0 (Eq. S9 from Evers et al., 2010 JxBot,61:2203–2216)

    VPD = utils.vapor_pressure_deficit(Tac, Tlc, hs)

    x1c,x2c,x1j,x2j,x1t,x2t,Rd = compute_an_2par(par_photo, PPFD, Tlc)

    An, Cc, Ci, gs = compute_an_analytic(Tlc,VPD,x1c,x2c,x1j,x2j,x1t,x2t,Rd,psi,model,g0,rbt,ca,m0,psi0,D0,n)

    return An, Cc, Ci, gs


def transpiration_rate(Tlc, ea, gs, gb, Pa = 101.3):
    """
    Returns E, leaf transpiration in [mol m-2 s-1]
    
    :Parameters:
    - **Tlc**: float, leaf temperature [degreeC]
    - **ea**: float, air vapor pressure [kPa]
    - **gs**: stomatal conductance for water [mol m-2 s-1]
    - **gb**: float, boundary layer conductance for H2O [mol m-2 s-1]
    - **Pa**: float, atmospheric pressure [kPa].
    
    :Worning:
    gb is calculated for both sides of a flat leaf.
    """

#    gva = gb*1.37 # boundary layer conductance for water vapour transport [mol m2 s-1] # R: ancienne formule : gb*1.4
    gv = 1./((2./gb)+(1./gs)) # Formulation by Pearcy 1989
    es_l = utils.saturated_air_vapor_pressure(Tlc)  # Saturated vapor pressure at leaf surface [kPa]
    E = gv*((es_l-ea)/Pa) # Transpiration rate [mol m-2 s-1]

    return E


def gas_exchange_rates(g, par_photo, par_photo_N, par_gs, meteo, E_type2,
                       leaf_lbl_prefix='L', rbt=2./3.):
    """
    Calculates gas exchange fluxes at the leaf scale according to the analytical scheme described by Evers et al. (JxBot 2010, 2203–2216).

    :Parameters:
    - **g**: a multiscale tree graph object
    - **par_photo**: dictionary, the parameters of the Farquhar model for net CO2 assimilation (cd :func:`par_photo_default`)
    - **par_gs**: dictionary, the parameters of the stomatal conductance model (model, g0, m0, psi0, D0, n)
    - **meteo**: a `pandas.DataFrame`-like object
    - **E_type2**: string, one of two 'Ei' or 'Eabs'
    - **leaf_lbl_prefix**: string, the prefix of the label of the leaves
    - **rbt**: float, the combined turbulance and boundary layer resistance for CO2 transport [m2 s ubar umol-1]
    - **ca**: float, CO2 partial pressure of the air [ubar]
    
    
    :Attaches to each leaf:
    - **An**: the net CO2 assimilation [umol m-2 s-1]
    - **Ci**: intercellular CO2 concentration [umol mol]
    - **gs**: stomatal conductance to water vapor [mol m-2 s-1]
    - **gb**: boundary layer conductance to water vapor [mol m-2 s-1]
    - **E** leaf transpiration [mol m-2 s-1].
    """
    
    model,g0max,m0,psi0,D0,n = [par_gs[ikey] for ikey in ('model','g0','m0','psi0','D0','n')]
    
    meteo_leaf = deepcopy(meteo)
    meteo_leaf = meteo_leaf.iloc[0]
    
    for vid in g:
        if vid>0:
            node=g.node(vid)
            if node.label.startswith(leaf_lbl_prefix):
                Tac = meteo_leaf.Tac
                hs = meteo_leaf.hs
                u = meteo_leaf.u
                Ca = meteo_leaf.Ca
                Pa = meteo_leaf.Pa
                
                node.u = u  # TODO replace the meso-wind speed (u) by a micro-wind speed at the level of each leaf

                psi = node.properties()['psi_head'] #leaf water potential [MPa] (assumed equal to that of the petiole)
                Tlc = node.properties()['Tlc']
                PPFD_leaf = node.properties()[E_type2]
#                w = node.Length / 100. # leaf length in m                

                meteo_leaf['PPFD'] = PPFD_leaf
                meteo_leaf['Rg'] = PPFD_leaf/(0.48*4.6)
#                meteo_leaf['u'] = min(1., meteo_leaf['u']) # Hack: see Nobel p.338
                
#                if not 'par_photo' in node.properties():
                leaf_par_photo = deepcopy(par_photo)
                leaf_par_photo['Vcm25'] = par_photo_N['Vcm25_N'][0]*node.Na+par_photo_N['Vcm25_N'][1]
                leaf_par_photo['Jm25'] = par_photo_N['Jm25_N'][0]*node.Na+par_photo_N['Jm25_N'][1]
                leaf_par_photo['TPU25'] = par_photo_N['TPU25_N'][0]*node.Na+par_photo_N['TPU25_N'][1]
                leaf_par_photo['Rd'] = par_photo_N['Rd_N'][0]*node.Na+par_photo_N['Rd_N'][1]
                dHd_max = leaf_par_photo['dHd']
                dHd = dHd_sensibility(psi, Tlc, dhd_max=dHd_max, dhd_inhib_beg=195., dHd_inhib_max=190.,
                                      psi_inhib_beg=-.75, psi_inhib_max=-2., temp_inhib_beg=32, temp_inhib_max=33)

                leaf_par_photo['dHd'] = dHd
                node.par_photo = leaf_par_photo
                
                g0 = g0max#*g0_sensibility(psi, psi_crit=-1, n=4)

                An, Cc, Ci, gs = an_gs_ci(node.par_photo, meteo_leaf, psi, Tlc,
                                          model, g0, rbt, Ca, m0, psi0, D0, n)

                # Boundary layer conductance
                # l_w = node.Length/100.*0.72 # leaf length in the downwind direction [m]
                # d_bl = 4.*(l_w/max(1.e-3,u))**0.5 /1000. # Boundary layer thikness in [m] (Nobel, 2009 pp.337)
                # Dj0 = 2.13*1.e-5 #[m2 s-1] at P=1. atm and t=0. °C (Nobel, pp.545)
                # Dj = Dj0*(101.3/Pa)*((Tac+273.15)/273.15)**1.8 #(Nobel, eq.8.8, pp.379)
                # gb = Dj*(Pa * 1.e-3) / ((R*1.e-3) * (Tac+273.15) * d_bl) # [mol m-2 s-1] (Nobel, 2009 pp.337)
                gb = boundary_layer_conductance(node.Length, u, Pa, Tac, R)

                # Transpiration
                es_a = utils.saturated_air_vapor_pressure(Tac)
                ea = es_a*hs/100.
                E = transpiration_rate(Tlc, ea, gs, gb, Pa)

                node.An = An
                node.Ci = Ci
                node.gs = gs
                node.gb = gb
                node.E = max(0.,E)

    return