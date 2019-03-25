"""
@author: Gaetan Louarn
@author: Rami ALBASHA

This is a library for calculating meteorological variables
"""

from math import *
from scipy import exp

R = 8.314510  # L kPa mol-1 K-1


def saturated_air_vapor_pressure(Tac):
    """Compute saturated air vapor pressure.

    Args:
        Tac (float): [degreeC] air temperature

    Returns:
        (float): [kPa] saturated air vapor pressure

    """
    return 0.611 * exp(17.27 * Tac / (237.3 + Tac))


def relative_humidity(ea, es_a):
    """Compute air relative humidity

    Args:
        ea (float): [kPa] actual air vapor pressure
        es_a (float): [kPa] saturated air vapor pressure

    Returns:
        (float): [-] air relative humidity (%, between 0 and 100)

    """
    return (ea / es_a) * 100.


def air_vapor_pressure_deficit(Tac, hs):
    """Compute air vapore pressure deficit

    Args:
        Tac (float): [degreeC] air temperature
        hs (float): (%) air relative humidity

    Returns:
        (float): [kPa] air vapor pressure deficit

    """
    es_a = saturated_air_vapor_pressure(Tac)
    ea = es_a * hs / 100
    return es_a - ea


def solar_declination(DOY):
    """Compute declination angle of the sun

    Args:
        DOY (float): [-] day of year

    Returns:
        (float): [rad] solar declination

    """
    alpha = 2 * 3.14 * (DOY - 1) / 365
    return (0.006918 - 0.399912 * cos(alpha) + 0.070257 * sin(alpha))


def extraterrestrial_solar_irradiance(DOY, HU, latitude):
    """ Compute extraterrestrial irradiance flux density

    Args:
        DOY (float): [-] day of year
        HU (float): [h] hour of the day
        latitude (float): [degree] location latitude

    Returns:
        (float): [W m-2] solar extraterrestrial irradiance flux density

    """
    hrad = 2 * 3.14 / 24 * (HU - 12)
    lat = radians(latitude)
    dec = solar_declination(DOY)
    costheta = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrad)
    Io = 1370 * (1 + 0.033 * cos(2 * 3.14 * (DOY - 4) / 366))
    return Io * costheta
    # extraterrestrial_solar_irradiance (100,11,44)


def clear_sky_global_radiation(Ra, z):
    """ compute clear sky global radiation according to extraterestrial radiation and altitude (m)
    eq 37 - FA056, p 51"""
    return Ra * (3600. / 1e6) * (0.75 + 2e-5 * z)  # formule pour Ra en MJ.m-2.h-1
    # clear_sky_global_radiation(extraterrestrial_solar_irradiance (100,11,44),10.)


def canopy_net_longwave_radiation_loss(Tac, Rs, Rs0):
    """ compute LongwaveNetRadiation = f(Tk_a,Tac_x, ea, Rs,Rs0) after FA0 56 eq #(39)
    Tak  : absolute air temperature (K)
    es_a : saturated vapor pressure in the ambiant air (kPa)
    Rs : measured solar radiation (MJ m2 hour)
    Rs0 : calculated clear-sky solar radiation (MJ m2 hour)"""

    ## Declaration des constantes
    sigma = 2.0412 * 1e-10  # Stefan-Boltzmann constant per surface area (MJ m-2 K-4 H-1)

    ## Compute longwave net radiation
    es_a = saturated_air_vapor_pressure(Tac)  # % saturated vapor pressure in the ambiant air (kPa)
    Tak = Tac + 273.16
    ratio = min(1., Rs / Rs0)
    Rnl = sigma * (Tak ** 4) * (0.34 - 0.14 * (es_a ** 0.5)) * (1.35 * ratio - 0.35)
    return Rnl * 1e6 / 3600.  # (en w.m-2)
    # canopy_net_longwave_radiation_loss(25., 600.*3600./1e6, clear_sky_global_radiation(extraterrestrial_solar_irradiance(100,11,44),10.))


def net_absorbed_radiation(Rg, Tac, DOY, HU, latitude=0.44, altitude=0., albedo=0.2):
    """ """
    Rextra = extraterrestrial_solar_irradiance(DOY, HU, latitude)
    Rs0_ = clear_sky_global_radiation(Rextra, altitude)
    Rnl = canopy_net_longwave_radiation_loss(Tac, Rg * 3600. / 1e6, Rs0_)
    Rabs = Rg * (1 - albedo) - Rnl

    return Rabs
    # net_absorbed_radiation(600., 25., 100, 11)


def boundary_layer_conductance(u, w=0.1):
    """ compute Boundary layer conductance for CO2= f(u, w) after Kim and Lieth (2003)
    u : wind speed m s-1
    w : leaf Characteristic dimension in relation to wind speed (m)"""

    d = 0.72 * w  # leaf dimension (m)
    gb = 0.147 * (u / d) ** 0.5  # boundary layer conductance (mol m2 s-1)
    # gb = 3.33 #  pour les simulations a niveau de feuille, il faut introduire la valeur de gb du LcPro (rb =0.33 mol m-2 s-1)
    return gb
    # boundary_layer_conductance(2.)


def celsius_to_kelvin(T):
    """
    Converts from Celsius to absolute temperature.
    """
    Tak = T + 273.
    return Tak


def vapor_pressure_deficit(Tac, Tlc, hs):
    """
    Returns leaf to air vapour pressure deficit [kPa].
    
    Parameters:
    -`Tac`: air temperature [degreeC]
    -`Tlc`: leaf temperature [degreeC]
    - `hs`: air relative humidity (%)
    """
    es_l = saturated_air_vapor_pressure(Tlc)  # % saturated vapor pressure in the leaf (kPa)
    es_a = saturated_air_vapor_pressure(Tac)  # % saturated vapor pressure in the ambiant air (kPa)
    ea = es_a * hs / 100  # % vapor pressure in the ambiant air (kPa)
    return es_l - ea


def cmol2cpa(temp, conc=400.):
    """
    Returns CO2 partial pressure [ubar] as a function of temperature [degrees]
    """
    V = 1.e6 * R * (temp + 273) / 101.3  # Volume of 1000000 mols of air [L]
    P = conc * R * (temp + 273) / V  # CO2 partial pressure [kPa]
    return P * 1.e4


def cpa2cmol(temp, partpress):
    """
    Returns CO2 concentration [umol mol-1] as a function of temperature [degrees]
    """
    V = 1.e6 * R * (temp + 273) / 101.3  # Volume of 1000000 mols of air [L]
    conc = partpress * V / (R * (temp + 273))  # CO2 partial pressure [kPa]
    return conc * 1.e-4
