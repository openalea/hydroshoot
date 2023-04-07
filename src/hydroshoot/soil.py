from math import ceil, exp, pi, log
from typing import List

from numpy import mean

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


def calc_volumetric_water_content_from_water_reservoir(water_content: float, soil_thickness: float) -> float:
    """Computes soil volumetric water content

    Args:
        water_content: [mm] soil water content
        soil_thickness: [mm] soil layer thickness

    Returns:
        (float): [-] soil volumetric water content

    """

    return water_content / soil_thickness


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
    theta = max(theta, theta_res * (1 + 1.e-6))
    m = 1 - 1. / n
    if theta == theta_sat:
        psi_soil = 0
    else:
        s_e = (theta - theta_res) / (theta_sat - theta_res)
        psi_soil = - 1. / alpha * ((1. / s_e) ** (1. / m) - 1) ** (1. / n)

    return psi_soil


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
    effective_saturation = 1. / ((1. + abs(alpha * psi) ** n) ** m)
    return k_sat * (effective_saturation ** 0.5) * (1. - (1. - effective_saturation ** (1. / m)) ** m) ** 2


def calc_root_soil_resistance(soil_conductivity: float, root_radius: float = 0.0001,
                              root_length_density: float = 2000) -> float:
    """Calculates the resistance to water flow at the soil-root interface according to Gardner (1960).

    Args:
        root_radius: [m] average radius of roots
        soil_conductivity: [cm d-1] soil hydraulic conductivity
        root_length_density: [m m-3] root length per soil volume

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
    k = soil_conductivity * 10. / 86400.  # cm d-1 -> kg m-2 s-1
    d = (pi * root_length_density) ** -0.5  # half distance between roots [m]
    return log(d ** 2 / root_radius ** 2) / (4 * pi * k)  # m2 s kg-1


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
        root_radius=root_radius,
        root_length_density=root_length / rhyzosphere_volume)
    uptake_per_unit_root_length = transpiration / root_length  # kg m-1 s-1
    return bulk_soil_water_potential - resistance * uptake_per_unit_root_length * (
            cst.water_density * cst.gravitational_acceleration * 1.e-6)  # m -> MPa


class SoilTexture:
    def __init__(self, theta_res: float, theta_sat: float, alpha: float, n: float, k_sat: float):
        self.theta_res = theta_res
        self.theta_sat = theta_sat
        self.alpha = alpha
        self.n = n
        self.k_sat = k_sat

    def to_dict(self):
        return {k: getattr(self, k) for k in ('theta_res', 'theta_sat', 'alpha', 'n', 'k_sat')}


class SoilVoxel:
    def __init__(self, theta_init: float, thickness: float, theta_res: float, theta_fc: float, theta_sat: float,
                 alpha: float, n: float):
        self.theta = theta_init
        self.thickness = thickness
        self.theta_res = theta_res
        self.theta_fc = theta_fc
        self.theta_sat = theta_sat
        self.alpha = alpha
        self.n = n

        self.water_content_fc = self.theta_fc * self.thickness
        self.water_content_res = self.theta_res * self.thickness
        self.water_content = self.theta * self.thickness
        self.water_potential = calc_soil_water_potential(
            theta=self.theta,
            theta_res=self.theta_res,
            theta_sat=self.theta_sat,
            alpha=self.alpha,
            n=self.n)

        self.drainage = None
        self.root_density = None

    def calc_water_balance(self, water_input: float, water_uptake: float):
        water_content = max(self.water_content_res, self.water_content + water_input - water_uptake)
        self.drainage = max(0.0, water_content - self.water_content_fc)
        self.water_content = water_content - self.drainage
        self.theta = calc_volumetric_water_content_from_water_reservoir(
            water_content=self.water_content,
            soil_thickness=self.thickness)
        self.water_potential = calc_soil_water_potential(
            theta=self.theta,
            theta_res=self.theta_res,
            theta_sat=self.theta_sat,
            alpha=self.alpha,
            n=self.n)


class Soil:
    def __init__(self, texture_class_name: str, depth: float, thickness_single_layer: float,
                 root_density_model: str, root_density_decrease_shape_param: float,
                 root_density_at_max_depth: float, water_potential_init: float = None):

        psi_fc = -0.033
        self.water_potential = min(psi_fc, water_potential_init) if water_potential_init is not None else psi_fc
        self.texture_class_name = texture_class_name
        self.depth = depth
        self.voxels = self.set_voxels(
            thickness_single_layer=thickness_single_layer,
            water_potential_init=self.water_potential)
        self.set_root_density_profile(
            root_density_model=root_density_model,
            root_density_decrease_shape_param=root_density_decrease_shape_param,
            root_density_at_max_depth=root_density_at_max_depth)

    def set_voxels(self, thickness_single_layer: float, water_potential_init: float):
        nb_layers = ceil(self.depth / thickness_single_layer)

        soil_class = SoilTexture(*SOIL_PROPS[self.texture_class_name])
        theta_fc = calc_volumetric_water_content_from_water_potential(
            psi=-33 * 10,  # kPa to cm
            theta_res=soil_class.theta_res,
            theta_sat=soil_class.theta_sat,
            alpha=soil_class.alpha,
            n=soil_class.n)
        theta_init = calc_volumetric_water_content_from_water_potential(
            psi=water_potential_init * 1.e4,  # MPa to cm
            theta_res=soil_class.theta_res,
            theta_sat=soil_class.theta_sat,
            alpha=soil_class.alpha,
            n=soil_class.n)

        voxels = []
        depth_upper = 0
        for i in range(nb_layers):
            voxels.append(SoilVoxel(
                theta_init=theta_init,
                thickness=min(thickness_single_layer, self.depth - depth_upper),
                theta_res=soil_class.theta_res,
                theta_fc=theta_fc,
                theta_sat=soil_class.theta_sat,
                alpha=soil_class.alpha,
                n=soil_class.n))

        return voxels

    def set_root_density_profile(self, root_density_model: str, root_density_decrease_shape_param: float,
                                 root_density_at_max_depth: float):
        root_density_profile = calc_root_density_profile(
            model=root_density_model,
            a=root_density_decrease_shape_param,
            shift=root_density_at_max_depth,
            soil_layers=[vox.thickness for vox in self.voxels])
        for vox, root_density in zip(self.voxels, root_density_profile):
            vox.root_density = root_density

    def calc_water_balance(self, water_input: float, transpiration: float):
        water_input_profile = [water_input] + [0] * (len(self.voxels) - 1)
        water_uptake_profile = [transpiration * vox.root_density for vox in self.voxels]
        drainage = 0
        for i, vox in enumerate(self.voxels):
            vox.calc_water_balance(
                water_input=water_input_profile[i] + drainage,
                water_uptake=water_uptake_profile[i])
            drainage = vox.drainage

    def calc_lumped_water_potential(self):
        self.water_potential = mean([vox.root_density * vox.water_potential for vox in self.voxels])


def calc_root_density_profile(model: str, a: float, shift: float, soil_layers: List[float]) -> List[float]:
    soil_depth = sum(soil_layers)

    depth_upper = 0
    root_density = []
    integral = lambda x: ((a / 2 * (x ** 2) if model == 'linear' else (1 / a * exp(-a * (soil_depth - x)))) + shift * x)
    for layer_thickness in soil_layers:
        root_density.append(
            integral(soil_depth - depth_upper) - integral(soil_depth - (depth_upper + layer_thickness)))
        depth_upper += layer_thickness

    density_tot = sum(root_density)

    return [v / density_tot for v in root_density]
