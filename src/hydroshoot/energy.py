# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Energy balance module of HydroShoot.

This module computes leaf (and eventually other elements) tempertaure of a
given plant shoot.
"""

from math import pi

import alinea.astk.icosphere as ico
import openalea.plantgl.all as pgl
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import turtle
from openalea.mtg.mtg import MTG
from scipy import optimize
from sympy import Symbol
from sympy.solvers import nsolve

import hydroshoot.constants as cst
from hydroshoot import utilities as utils
from hydroshoot.architecture import get_leaves


def pgl_scene(g, flip=False):
    geometry = g.property('geometry')
    scene = pgl.Scene()
    for id in geometry:
        if not flip:
            sh = pgl.Shape(geometry[id])
        else:
            sh = pgl.Shape(pgl.AxisRotated(pgl.Vector3(1, 0, 0), pi, geometry[id]))
        sh.id = id
        scene.add(sh)
    return scene


def set_form_factors_simplified(g, pattern=None, infinite=False, leaf_lbl_prefix='L', turtle_sectors='46',
                                icosphere_level=3, unit_scene_length='cm'):
    """Computes sky and soil contribution factors (resp. k_sky and k_soil) to the energy budget equation.
    Both factors are calculated and attributed to each element of the scene.

    Args:
        g: a multiscale tree graph object
        pattern (tuple): 2D Coordinates of the domain bounding the scene for its replication.
            (xmin, ymin, xmax, ymax) scene is not bounded along z axis.
            Alternatively a *.8 file.
            if `None` (default), scene is not repeated
        infinite (bool): Whether the scene should be considered as infinite
            (see :func:`runCaribu` from `CaribuScene` package)
        leaf_lbl_prefix (str): the prefix of the leaf label
        turtle_sectors (str): number of turtle sectors
            (see :func:`turtle` from `sky_tools` package)
        icosphere_level (int): the level of refinement of the dual icosphere
            (see :func:`alinea.astk.icosphere.turtle_dome` for details)
        unit_scene_length (str): the unit of length used for scene coordinate and for pattern
            (should be one of `CaribuScene.units` default)

    Returns:
        g updated with the properties 'ff-sky', 'ff_leaves' and 'ff_soil' (resp. for sky, leaves and soil form factors)

    Notes:
        This function is a simplified approximation of the form factors matrix which is calculated by the
            function :func:`form_factors_matrix`. The canopy is turned upside-down and light is projected in each
            case to estimate the respective contribution of the sky ({z}>=0) and soil ({z}<=0) to energy budget
            calculations. This is printed as pirouette-cacahuete in the console!
        When **icosphere_level** is defined, **turtle_sectors** is ignored.

    """
    geom = g.property('geometry')
    label = g.property('label')
    opts = {'SW': {vid: ((0.001, 0) if label[vid].startswith(leaf_lbl_prefix) else (0.001,)) for vid in geom}}
    if not icosphere_level:
        direction = turtle.turtle(sectors=turtle_sectors, format='uoc', energy=1.)[2]
    else:
        vert, fac = ico.turtle_dome(icosphere_level)
        direction = ico.sample_faces(vert, fac, iter=None, spheric=False).values()
        direction = [i[0] for i in direction]
        direction = map(lambda x: tuple(list(x[:2]) + [-x[2]]), direction)

    caribu_source = list(zip(len(direction) * [1. / len(direction)], direction))
    k_soil, k_sky, k_leaves = {}, {}, {}

    for s in ('pirouette', 'cacahuete'):
        print('... %s' % s)
        if s == 'pirouette':
            scene = pgl_scene(g, flip=True)
        else:
            scene = pgl_scene(g)

        caribu_scene = CaribuScene(scene, light=caribu_source, opt=opts, scene_unit=unit_scene_length, pattern=pattern)

        # Run caribu
        raw, aggregated = caribu_scene.run(direct=True, infinite=infinite, split_face=False, simplify=True)

        if s == 'pirouette':
            k_soil_dict = aggregated['Ei']
            max_k_soil = float(max(k_soil_dict.values()))
            k_soil = {vid: k_soil_dict[vid] / max_k_soil for vid in k_soil_dict}
        elif s == 'cacahuete':
            k_sky_dict = aggregated['Ei']
            max_k_sky = float(max(k_sky_dict.values()))
            k_sky = {vid: k_sky_dict[vid] / max_k_sky for vid in k_sky_dict}

    for vid in aggregated['Ei']:
        k_leaves[vid] = max(0., 2. - (k_soil[vid] + k_sky[vid]))

    g.properties()['ff_sky'] = k_sky
    g.properties()['ff_leaves'] = k_leaves
    g.properties()['ff_soil'] = k_soil
    return g


def set_leaf_temperature_to_air_temperature(g: MTG, leaf_lbl_prefix: str):
    """Basic model for leaf temperature, considered equal to air temperature for all leaves

    Args:
        g: a multiscale tree graph object
        leaf_lbl_prefix (str): the prefix of the leaf label

    Returns:
        the mtg updated with all leaves having their temperature equal to that of the (local) air

    """
    g.properties()['Tlc'] = {vid: g.node(vid).Tac for vid in get_leaves(g=g, leaf_lbl_prefix=leaf_lbl_prefix)}
    return g


def set_local_wind_speed(g, meteo, leaf_lbl_prefix='L') -> dict:
    """Basic model for wind speed at leaf level, considered equal to air wind speed for all leaves

    Args:
        g: a multiscale tree graph object
        meteo (DataFrame): forcing meteorological variables
        leaf_lbl_prefix (str): the prefix of the leaf label

    Returns:
        (dict): keys are leaves vertices ids and their values are all equal to air wind speed

    """
    leaves = get_leaves(g, leaf_lbl_prefix)
    u = meteo.u[0]
    return {vid: u for vid in leaves}


def set_local_air_temperature(g, meteo, leaf_lbl_prefix='L') -> dict:
    """Basic model for air temperature at the leaf level, considered equal to air temperature at reference height for
    all leaves.

    Args:
        g: a multiscale tree graph object
        meteo (DataFrame): forcing meteorological variables
        leaf_lbl_prefix (str): the prefix of the leaf label

    Returns:
        (dict): keys are leaves vertices ids and their values are all equal to air temperature at reference height

    """
    return {vid: meteo['Tac'].iloc[0] for vid in get_leaves(g=g, leaf_lbl_prefix=leaf_lbl_prefix)}


def set_local_vpd(g: MTG, relative_humidity: float, leaf_lbl_prefix: str) -> MTG:
    """Calculates leaf-to-air vapor pressure deficit for all leaves.

    Args:
        g: a multiscale tree graph object
        relative_humidity: (%) air relative humidity (between 0 and 100)
        leaf_lbl_prefix: the prefix of the leaf label

    Returns:
        the mtg updated with vapor pressure deficit values

    """
    g.properties()['vpd'] = {vid: utils.vapor_pressure_deficit(temp_air=g.node(vid).Tac,
                                                               temp_leaf=g.node(vid).Tlc,
                                                               rh=relative_humidity)
                             for vid in get_leaves(g=g, leaf_lbl_prefix=leaf_lbl_prefix)}
    return g


def _gbH(leaf_length, wind_speed):
    """Computes boundary layer conductance to heat

    Args:
        leaf_length (float): [m] leaf length
        wind_speed (float): [m s-1] local wind speed

    Returns:
        (float): [W m-2 K-1] boundary layer conductance to heat

    References:
        Nobel P. 2005.
            Temperature and energy budgets.
            In Nobel S, eds. Physicochemical and Environmental Plant Physiology.
            Elsevier Academic Press, 307–350.

    """
    l_w = leaf_length * 0.72  # leaf length in the downwind direction [m]
    d_bl = 4. * (l_w / max(1.e-3, wind_speed)) ** 0.5 / 1000.  # Boundary layer thickness in [m] (Nobel, 2009 pp.337)
    return 2. * 0.026 / d_bl  # Boundary layer conductance to heat [W m-2 K-1]


def calc_heat_boundary_layer_conductance(g, leaf_ids, unit_scene_length):
    """Calculates boundary conductance to heat transfer for each leaf.

    Args:
        g: mtg object
        leaf_ids (list of int): leaf ids
        unit_scene_length (str): the unit of length used for scene coordinate and for pattern
            (should be one of `CaribuScene.units` default)

    Returns:
        g with the boundary layer conductance to heat 'gbH' in [W m-2 K-1] property for each leaf

    """
    conv = {'mm': 1.e-3, 'cm': 1.e-2, 'm': 1.}[unit_scene_length]
    for vid in leaf_ids:
        g.node(vid).gbH = _gbH(
            leaf_length=g.node(vid).properties()['Length'] * conv,
            wind_speed=g.node(vid).properties()['u'])
    return g


def calc_leaf_temperature(g, t_soil, t_sky_eff, leaf_ids, max_iter=100, t_error_crit=0.01, t_step=0.5):
    """Computes the temperature of each individual leaf and soil elements.

    Args:
        g: a multiscale tree graph object
        t_soil (float): [°C] soil surface temperature
        t_sky_eff (float): [°C] effective sky temperature
        solo (bool):
            if True (default), calculates energy budget for each element assuming the temperatures of surrounding
                leaves as constant (from previous calculation step)
            if False, computes simultaneously all temperatures using `sympy.solvers.nsolve` (**very costly!!!**)
        leaf_ids (list of int): leaf ids
        max_iter (int): maximum allowed iteration (used only when :arg:`solo` is True)
        t_error_crit (float): [°C] maximum allowed error in leaf temperature (used only when :arg:`solo` is True)
        t_step (float): [°C] maximum temperature step between two consecutive iterations

    Returns:
        (dict): [°C] the tempearture of individual leaves given as the dictionary keys
        (int): [-] the number of iterations (not None only when :arg:`solo` is True)

    """

    temp_sky = utils.celsius_to_kelvin(t_sky_eff)
    temp_soil = utils.celsius_to_kelvin(t_soil)

    t_prev = g.property('Tlc')
    t_new = None

    t_error_trace = []
    it = 0
    for it in range(max_iter):
        t_new = {}

        for vid in leaf_ids:
            mtg_node = g.node(vid)
            shortwave_inc = mtg_node.properties()['Rg']
            ff_sky = mtg_node.properties()['ff_sky']
            ff_leaves = mtg_node.properties()['ff_leaves']
            ff_soil = mtg_node.properties()['ff_soil']
            gb_h = mtg_node.properties()['gbH']
            evap = mtg_node.properties()['E']
            t_leaf = utils.celsius_to_kelvin(mtg_node.properties()['Tlc'])
            temp_air = utils.celsius_to_kelvin(mtg_node.properties()['Tac'])

            longwave_gain_from_leaves = ff_leaves * cst.sigma * t_leaf ** 4

            def _VineEnergyX(t_leaf):
                shortwave_abs = cst.a_glob * shortwave_inc
                longwave_net = (cst.e_leaf * (ff_sky * cst.e_sky * cst.sigma * temp_sky ** 4 +
                                              cst.e_leaf * longwave_gain_from_leaves +
                                              ff_soil * cst.e_soil * cst.sigma * temp_soil ** 4)
                                - 2 * cst.e_leaf * cst.sigma * t_leaf ** 4)
                latent_heat_loss = -cst.lambda_ * evap
                sensible_heat_net = -gb_h * (t_leaf - temp_air)
                energy_balance = shortwave_abs + longwave_net + latent_heat_loss + sensible_heat_net
                return energy_balance

            t_new[vid] = utils.kelvin_to_celsius(optimize.newton_krylov(_VineEnergyX, t_leaf))

        t_error = max([abs(t_prev[vtx] - t_new[vtx]) for vtx in leaf_ids])
        t_error_trace.append(t_error)

        if t_error < t_error_crit:
            break
        else:
            try:
                if abs(t_error_trace[-1] - t_error_trace[-2]) < t_error_crit:
                    t_step = max(0.01, t_step / 2.)
            except IndexError:
                pass

            t_prev = {k: t_prev[k] + t_step * (t_new[k] - t_prev[k]) for k in leaf_ids}

    return t_new, it


def calc_leaf_temperature2(g, t_soil, t_sky_eff, leaf_ids):
    """Computes the temperature of each individual leaf and soil elements whilst considering explicitly each other
    leaf temperature energy in a matrix solution using symbolic solver.

    Args:
        g: a multiscale tree graph object
        t_soil (float): [°C] soil surface temperature
        t_sky_eff (float): [°C] effective sky temperature
        leaf_ids (list of int): leaf ids

    Returns:
        (dict): [°C] the temperature of individual leaves given as the dictionary keys
        (int): [-] the number of iterations (not None only when :arg:`solo` is True)

    """

    temp_sky = utils.celsius_to_kelvin(t_sky_eff)
    temp_soil = utils.celsius_to_kelvin(t_soil)

    t_prev = g.property('Tlc')

    t_lst = []
    t_new = {vid: Symbol('t%d' % vid) for vid in leaf_ids}

    eq_lst = []
    t_leaf_lst = []
    for vid in leaf_ids:
        mtg_node = g.node(vid)
        shortwave_inc = mtg_node.properties()['Rg']

        ff_sky = mtg_node.properties()['ff_sky']
        ff_leaves = mtg_node.properties()['ff_leaves']
        ff_soil = mtg_node.properties()['ff_soil']
        temp_air = utils.celsius_to_kelvin(mtg_node.properties()['Tac'])

        gb_h = mtg_node.properties()['gbH']
        evap = mtg_node.properties()['E']
        t_leaf = t_prev[vid]

        t_leaf_lst.append(t_leaf)
        t_lst.append(t_new[vid])

        eq_aux = 0.
        for ivid in ff_leaves:
            if not g.node(ivid).label.startswith('soil'):
                eq_aux += -ff_leaves[ivid] * ((t_new[ivid]) ** 4)

        eq = (cst.a_glob * shortwave_inc +
              cst.e_leaf * cst.sigma * (ff_sky * cst.e_sky * (temp_sky ** 4) +
                                        cst.e_leaf * eq_aux + ff_soil * cst.e_soil * (temp_soil ** 4) -
                                        2 * (t_new[vid]) ** 4) -
              cst.lambda_ * evap - gb_h * (t_new[vid] - temp_air))

        eq_lst.append(eq)

    # tt = time.time()
    t_leaf0_lst = nsolve(eq_lst, t_lst, t_leaf_lst, verify=False) - 273.15
    # print("---%s seconds ---" % (time.time() - tt))

    t_new = {}
    for ivid, vid in enumerate(leaf_ids):
        t_new[vid] = float(t_leaf0_lst[ivid])
        ivid += 1

    return t_new, 1  # oneshot calculation (no iteration)


def force_soil_temperature(meteo):
    """A very simple model of soil temperature

    Args:
        meteo (DataFrame): forcing meteorological variables

    Returns:
        (double): [°C] soil temperature

    """
    if 'Tsoil' in meteo.columns:
        t_soil = meteo.Tsoil[0]
    else:
        dt_soil = [3, 3, 3, 3, 3, 3, 3, 3, 10, 15, 20, 20, 20, 20, 20, 15, 6, 5, 4, 3, 3, 3, 3, 3]
        t_soil = meteo.Tac[0] + dt_soil[meteo.index.hour[0]]
    return t_soil


def calc_effective_sky_temperature(diffuse_to_total_irradiance_ratio: float, temperature_cloud: float,
                                   temperature_sky: float) -> float:
    # Change the t_sky_eff formula (cf. Gliah et al., 2011, Heat and Mass Transfer, DOI: 10.1007/s00231-011-0780-1)
    return diffuse_to_total_irradiance_ratio * temperature_cloud + (
            1 - diffuse_to_total_irradiance_ratio) * temperature_sky
