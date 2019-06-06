# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Some useful common functions.
"""

from math import exp

ideal_gas_cst = 8.314510  # L kPa mol-1 K-1
absolute_zero = -273.15  # absolute zero temperature


def saturated_air_vapor_pressure(temp):
    """Compute saturated air vapor pressure.

    Args:
        temp (float): [°C] air temperature

    Returns:
        (float): [kPa] saturated air vapor pressure

    """
    return 0.611 * exp(17.27 * temp / (237.3 + temp))


def celsius_to_kelvin(temp):
    """Converts temperature from Celsius to Kelvin.

    Args:
        temp (float): [°C] temperature

    Returns:
        (float): [K] temperature

    """

    return temp - absolute_zero


def kelvin_to_celsius(temp):
    """Converts temperature from Kelvin to Celsius.

    Args:
        temp (float): [K] temperature

    Returns:
        (float): [°C] temperature

    """

    return temp + absolute_zero


def vapor_pressure_deficit(temp_air, temp_leaf, rh):
    """Computes leaf-to-air vapour pressure deficit.

    Args:
        temp_air (float): [°C] air temperature
        temp_leaf (float): [°C] leaf temperature
        rh (float): [-] air relative humidity (%, between 0 and 1)

    Returns:
        (float): [kPa] leaf-to-air vapour pressure deficit

    """

    es_l = saturated_air_vapor_pressure(temp_leaf)  # % saturated vapor pressure in the leaf (kPa)
    es_a = saturated_air_vapor_pressure(temp_air)  # % saturated vapor pressure in the ambiant air (kPa)
    ea = es_a * rh / 100  # % vapor pressure in the ambiant air (kPa)

    return es_l - ea


def cmol2cpa(temp, concentration=400.):
    """Convert CO2 concentration into CO2 partial pressure

    Args:
        temp (float): [°C] leaf temperature
        partial_pressure (float): [ppm] CO2 concentration in the air

    Returns:
        (float): [ubar] CO2 partial pressure

    """

    volume = 1.e6 * ideal_gas_cst * celsius_to_kelvin(temp) / 101.3  # Volume of 1000000 mols of air [L]
    partial_pressure = concentration * ideal_gas_cst * celsius_to_kelvin(temp) / volume  # CO2 partial pressure [kPa]

    return partial_pressure * 1.e4


def cpa2cmol(temp, partial_pressure):
    """Convert CO2 partial pressure into CO2 concentration

    Args:
        temp (float): [°C] leaf temperature
        partial_pressure (float): [ubar] CO2 partial pressure

    Returns:
        (float): [umol mol-1] CO2 concentration

    """

    volume = 1.e6 * ideal_gas_cst * celsius_to_kelvin(temp) / 101.3  # Volume of 1000000 mols of air [L]
    conc = partial_pressure * volume / (ideal_gas_cst * celsius_to_kelvin(temp))  # CO2 partial pressure [kPa]

    return conc * 1.e-4
