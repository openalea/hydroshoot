from math import pi, log

from hydroshoot import constants as cst

SOIL_PROPS = dict(
    Sand=(0.045, 0.430, 0.145, 2.68, 712.8),
    Loamy_Sand=(0.057, 0.410, 0.124, 2.28, 350.2),
    Sandy_Loam=(0.065, 0.410, 0.075, 1.89, 106.1),
    Loam=(0.078, 0.430, 0.036, 1.56, 24.96),
    Silt=(0.034, 0.460, 0.016, 1.37, 6.00),
    Silty_Loam=(0.067, 0.450, 0.020, 1.41, 10.80),
    Sandy_Clay_Loam=(0.100, 0.390, 0.059, 1.48, 31.44),
    Clay_Loam=(0.095, 0.410, 0.019, 1.31, 6.24),
    Silty_Clay_Loam=(0.089, 0.430, 0.010, 1.23, 1.68),
    Sandy_Clay=(0.100, 0.380, 0.027, 1.23, 2.88),
    Silty_Clay=(0.070, 0.360, 0.005, 1.09, 0.48),
    Clay=(0.068, 0.380, 0.008, 1.09, 4.80))


def calc_volumetric_water_content_from_water_potential(
        psi: float, theta_res: float, theta_sat: float, alpha: float, n: float) -> float:
    """Computes soil water potential following van Genuchten (1980)

    Args:
        psi: [cm H2O] soil water potential
        theta_res: [m3(H2O) m-3(H2O)] soil residual volumetric water content
        theta_sat: [m3(H2O) m-3(H2O)] soil saturated volumetric water content
        alpha: [cm cm-1] shape parameter
        n: [-] shape parameter

    Returns:
        (float): [cm H2O] soil water potential

    Reference:
        van Genuchten M., 1980.
            A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
            Soil Science Society of America Journal 44, 892897.
    """

    m = 1 - 1 / n

    if psi == 0:
        theta = theta_sat
    else:
        theta = theta_res + (theta_sat - theta_res) / (1. + abs(alpha * psi) ** n) ** m

    return theta


def calc_soil_water_potential(theta: float, theta_res: float, theta_sat: float, alpha: float, n: float) -> float:
    """Computes soil water potential following van Genuchten (1980)

    Args:
        theta: [-] volumetric soil water content
        theta_res: [m3(H2O) m-3(H2O)] soil residual volumetric water content
        theta_sat: [m3(H2O) m-3(H2O)] soil saturated volumetric water content
        alpha: [cm cm-1] shape parameter
        n: [-] shape parameter

    Returns:
        (float): [cm H2O] soil water potential

    Reference:
        van Genuchten M., 1980.
            A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
            Soil Science Society of America Journal 44, 892897.
    """
    theta = min(max(theta, theta_res * (1 + 1.e-6)), theta_sat)
    m = 1 - 1. / n
    if theta == theta_sat:
        psi_soil = 0
    else:
        s_e = (theta - theta_res) / (theta_sat - theta_res)
        psi_soil = - 1. / alpha * ((1. / s_e) ** (1. / m) - 1) ** (1. / n)

    return psi_soil


def update_soil_water_potential(psi_soil_init, water_withdrawal, soil_class, soil_total_volume, psi_min):
    """Updates the value of soil water potential after subtracting water quantity taken by the plant.

    Args:
        psi_soil_init (float): [MPa] initial soil water potential
        water_withdrawal (float): [Kg T-1] water volume that is withdrawn from the soil (by transpiration for instance)
            during a timelapse T
        soil_class (str): one of the soil classes proposed by Carsel and Parrish (1988)
        soil_total_volume (float): [m3] total apparent volume of the soil
        psi_min (float): [MPa] minimum allowable water potential

    Returns:
        (float): [MPa] soil water potential calculated following the water retention curve of van Genuchten (1980)

    References:
        Carsel R., Parrish R., 1988.
            Developing joint probability distributions of soil water retention characteristics.
            Water Resources Research 24,755 - 769.
        van Genuchten M., 1980.
            A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
            Soil Science Society of America Journal 44, 892897.
    """

    psi_soil_init = min(-1e-6, psi_soil_init)
    theta_r, theta_s, alpha, n, k_sat = SOIL_PROPS[soil_class]

    m = 1. - 1. / n

    psi = psi_soil_init * 1.e6 / (cst.water_density * cst.gravitational_acceleration) * 100.  # MPa -> cm_H20
    theta_init = theta_r + (theta_s - theta_r) / (1. + abs(alpha * psi) ** n) ** m

    flux = water_withdrawal / cst.water_density  # kg T-1 -> m3 T-1

    porosity_volume = soil_total_volume * theta_s

    delta_theta = flux / porosity_volume  # [m3 m-3]

    theta = max(theta_r, theta_init - delta_theta)

    psi_soil = calc_soil_water_potential(theta=theta, theta_res=theta_r, theta_sat=theta_s, alpha=alpha, n=n) / (
            1.e6 / (cst.water_density * cst.gravitational_acceleration) * 100)

    return max(psi_min, psi_soil)


def calc_soil_conductivity(psi: float, soil_class: str) -> float:
    """Gives the actual soil hydraulic conductivity following van Genuchten (1980)

    Args:
        psi (float): [MPa] bulk soil water potential
        soil_class (str): soil texture classe according to Carsel and Parrish (1988)

    Returns:
        (float): [cm d-1] actual soil water conductivity

    References:
        Carsel R., Parrish R., 1988.
            Developing joint probability distributions of soil water retention characteristics.
            Water Resources Research 24,755 - 769.
        van Genuchten M., 1980.
            A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
            Soil Science Society of America Journal 44, 892897.
    """

    psi *= 1.e6 / (cst.water_density * cst.gravitational_acceleration) * 100.  # MPa -> cm_H20
    theta_r, theta_s, alpha, n, k_sat = SOIL_PROPS[soil_class]
    m = 1. - 1. / n
    effective_saturation = 1. / ((1. + (abs(alpha * psi)) ** n) ** m)
    return k_sat * (effective_saturation ** 0.5) * (1. - (1. - effective_saturation ** (1. / m)) ** m) ** 2


def calc_root_soil_resistance(soil_conductivity: float, rhyzosphere_volume: float,
                              root_radius: float = 0.0001, root_length: float = 2000) -> float:
    """Calculates the resistance to water flow at the soil-root interface according to Gardner (1960).

    Args:
        soil_conductivity: [cm d-1] soil hydraulic conductivity
        rhyzosphere_volume: [m3] volume of the rhyzosphere
        root_radius: [m] average radius of roots
        root_length: [m] root length

    Returns:
        (float): [m2 s kg-1] resistance to water flow at the soil-root interface

    References:
        Gardner (1960)
            Dynamic aspects of water availability to plants.
            Soil science 89, 63 - 73.
        Tardieu et al. (2015)
            Modelling the coordination of the controls of stomatal aperture, transpiration, leaf growth,
                and abscisic acid: update and extension of the Tardieu–Davies model.
            Journal of Experimental Botany 66, 2227 – 2237.
            https://doi.org/10.1093/jxb/erv039

    """
    root_length_density = root_length / rhyzosphere_volume  # [m m3]
    d = (pi * root_length_density) ** -0.5  # half distance between roots [m]
    k = soil_conductivity * 1.e-2 / 86400. * (2 * pi * d * root_length) * cst.water_density  # cm d-1 -> kg m-2 s-1
    return log(d ** 2 / ((2 * root_radius) ** 2)) / (4 * pi * k)  # m2 s kg-1


def calc_collar_water_potential(transpiration: float, bulk_soil_water_potential: float, rhyzosphere_volume: float,
                                soil_class: str, root_radius: float, root_length: float) -> float:
    """Calculates the lumped water potential of the root system, assumed to be equal to that at the plant collar.

    Args:
        rhyzosphere_volume: [m3] depth of the root system
        transpiration: [kg s-1] transpiration flux
        bulk_soil_water_potential: [MPa]
        soil_class (str): soil texture classe according to Carsel and Parrish (1988)
        root_radius: [m] average radius of roots
        root_length: [m] root length per plant

    Returns:
        [MPa] water potential at the plant collar

    """
    resistance = calc_root_soil_resistance(  # m2 s kg-1
        soil_conductivity=calc_soil_conductivity(
            psi=bulk_soil_water_potential,
            soil_class=soil_class),
        rhyzosphere_volume=rhyzosphere_volume,
        root_radius=root_radius,
        root_length=root_length)
    uptake_per_unit_root_length = transpiration / root_length  # kg m-1 s-1
    return bulk_soil_water_potential - resistance * uptake_per_unit_root_length * (
            cst.water_density * cst.gravitational_acceleration * 1.e-6)  # m -> MPa
