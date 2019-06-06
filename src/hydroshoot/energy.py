# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Energy balance module of HydroShoot.

This module computes leaf (and eventually other elements) tempertaure of a
given plant shoot.
"""

from scipy import optimize, mean
from sympy.solvers import nsolve
from sympy import Symbol
import time


from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import turtle
import alinea.astk.icosphere as ico
import openalea.plantgl.all as pgl
from math import pi

from hydroshoot import utilities as utils


def pgl_scene(g, flip=False):
    geometry = g.property('geometry')
    scene = pgl.Scene()
    for id in geometry:
        if not flip:
            sh = pgl.Shape(geometry[id])
        else:
            sh = pgl.Shape(pgl.AxisRotated(pgl.Vector3(1,0,0),pi,geometry[id]))
        sh.id = id
        scene.add(sh)
    return scene


def get_leaves(g, leaf_lbl_prefix='L'):
    label = g.property('label')
    return [vid for vid in g.VtxList() if
                             vid > 0 and label[vid].startswith(leaf_lbl_prefix)]


def get_leaves_length(g, leaf_lbl_prefix='L', length_lbl='Length', unit_scene_length='cm'):
    """get length of leaves of g [m]"""
    conv = {'mm': 1.e-3, 'cm': 1.e-2, 'm': 1.}[unit_scene_length]
    leaves = get_leaves(g, leaf_lbl_prefix)
    length = g.property(length_lbl)
    return {k: v * conv for k, v in length.iteritems() if k in leaves}



a_PAR = 0.87
a_NIR = 0.35
a_glob = 0.6
e_sky = 1.0
e_leaf = 0.96
e_soil = 0.95
sigma = 5.670373e-8
lambda_ = 44.0e3
Cp = 29.07


def form_factors_simplified(g, pattern=None, infinite=False, leaf_lbl_prefix='L', turtle_sectors='46',
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
        energy, emission, direction, elevation, azimuth = turtle.turtle(sectors=turtle_sectors, format='uoc', energy=1.)
    else:
        vert, fac = ico.turtle_dome(icosphere_level)
        direction = ico.sample_faces(vert, fac, iter=None, spheric=False).values()
        direction = [i[0] for i in direction]
        direction = map(lambda x: tuple(list(x[:2]) + [-x[2]]), direction)

    caribu_source = zip(len(direction) * [1. / len(direction)], direction)
    k_soil, k_sky, k_leaves = {}, {}, {}

    for s in ('pirouette', 'cacahuete'):
        print '... %s' % s
        if s == 'pirouette':
            scene = pgl_scene(g, flip=True)
        else:
            scene = pgl_scene(g)

        caribu_scene = CaribuScene(scene, light=caribu_source, opt=opts,
                                   scene_unit=unit_scene_length,
                                   pattern=pattern)

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

    return k_soil, k_sky, k_leaves


def leaf_temperature_as_air_temperature(g, meteo, leaf_lbl_prefix='L'):
    """Basic model for leaf temperature, considered equal to air temperature for all leaves

    Args:
        g: a multiscale tree graph object
        meteo (DataFrame): forcing meteorological variables
        leaf_lbl_prefix (str): the prefix of the leaf label

    Returns:
        (dict): keys are leaves vertices ids and their values are all equal to air temperature [°C]

    """
    leaves = get_leaves(g, leaf_lbl_prefix)
    t_air = meteo.Tac[0]
    return {vid: t_air for vid in leaves}


def leaf_wind_as_air_wind(g, meteo, leaf_lbl_prefix='L'):
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


def heat_boundary_layer_conductance(leaves_length, wind_speed=0):
    """Computes boundary layer conductance to heat of all mtg leaves

    Args:
        leaves_length (dict): [m] length of individual leaves given as the dictionary keys
        wind_speed (dict): [m s-1] local wind speed of individual leaves given as the dictionary keys

    Returns:
        (float): [W m-2 K-1] boundary layer conductance to heat of individual leaves given as the dictionary keys

    """

    if isinstance(wind_speed, dict):
        u = wind_speed
    else:
        u = {vid: wind_speed for vid in leaves_length}
    return {vid: _gbH(leaves_length[vid], u[vid]) for vid in leaves_length}


# TODO: split leaf_temperature() into two functions following whether solo is used or not
def leaf_temperature(g, meteo, t_soil, t_sky_eff, t_init=None, form_factors=None, gbh=None, ev=None, ei=None, solo=True,
                     ff_type=True, leaf_lbl_prefix='L', max_iter=100, t_error_crit=0.01, t_step=0.5):
    """Computes the temperature of each individual leaf and soil elements.

    Args:
        g: a multiscale tree graph object
        meteo (DataFrame): forcing meteorological variables
        t_soil (float): [°C] soil surface temperature
        t_sky_eff (float): [°C] effective sky temperature
        t_init(float or dict): [°C] temperature used for initialisation
        form_factors(3-tuple of float or dict): form factors for soil, sky and leaves.
            if None (default) (0.5, 0.5, 0.5) is used for all leaves
        gbh (float or dict): [W m-2 K-1] boundary layer conductance for heat
            if None (default) a default model is called with length=10cm and wind_speed as found in meteo
        ev (float or dict): [mol m-2 s-1] evaporation flux
            if None (default) evaporation is set to zero for all leaves
        ei (float or dict): [umol m-2 s-1] photosynthetically active radition (PAR) incident on leaves
            if None (default) PAR is set to zero for all leaves
        solo (bool):
            if True (default), calculates energy budget for each element assuming the temperatures of surrounding
                leaves as constant (from previous calculation step)
            if False, computes simultaneously all temperatures using `sympy.solvers.nsolve` (**very costly!!!**)
        ff_type (bool): form factor type flag. If true fform factor for a given leaf is expected to be a single value, or a dict of ff otherwxie
        leaf_lbl_prefix (str): the prefix of the leaf label
        max_iter (int): maximum allowed iteration (used only when :arg:`solo` is True)
        t_error_crit (float): [°C] maximum allowed error in leaf temperature (used only when :arg:`solo` is True)
        t_step (float): [°C] maximum temperature step between two consecutive iterations

    Returns:
        (dict): [°C] the tempearture of individual leaves given as the dictionary keys
        (int): [-] the number of iterations (not None only when :arg:`solo` is True)

    """

    leaves = get_leaves(g, leaf_lbl_prefix)
    it = 0

    if t_init is None:
        t_init = meteo.Tac[0]
    if form_factors is None:
        form_factors = 0.5, 0.5, 0.5
    if gbh is None:
        gbh = _gbH(0.1, meteo.u[0])
    if ev is None:
        ev = 0
    if ei is None:
        ei = 0

    k_soil, k_sky, k_leaves = form_factors
    properties = {}
    for what in ('t_init', 'gbh', 'ev', 'ei', 'k_soil', 'k_sky', 'k_leaves'):
        val = eval(what)
        if isinstance(val, dict):
            properties[what] = val
        else:
            properties[what] = {vid: val for vid in leaves}

    # macro-scale climatic data
    temp_sky = utils.celsius_to_kelvin(t_sky_eff)
    temp_air = utils.celsius_to_kelvin(meteo.Tac[0])
    temp_soil = utils.celsius_to_kelvin(t_soil)

    # initialisation
    t_prev = properties['t_init']

    # iterative calculation of leaves temperature
    if solo:
        t_error_trace = []
        it_step = t_step
        for it in range(max_iter):
            t_dict = {}

            for vid in leaves:
                shortwave_inc = properties['ei'][vid] / (0.48 * 4.6)  # Ei not Eabs

                ff_sky = properties['k_sky'][vid]
                ff_leaves = properties['k_leaves'][vid]
                ff_soil = properties['k_soil'][vid]

                gb_h = properties['gbh'][vid]
                evap = properties['ev'][vid]
                t_leaf = t_prev[vid]

                if not ff_type:
                    longwave_grain_from_leaves = -sigma * sum(
                        [ff_leaves[ivid] * (utils.celsius_to_kelvin(t_prev[ivid])) ** 4 for ivid in ff_leaves])
                else:
                    longwave_grain_from_leaves = ff_leaves * sigma * (utils.celsius_to_kelvin(t_leaf)) ** 4

                def _VineEnergyX(t_leaf):
                    shortwave_abs = a_glob * shortwave_inc
                    longwave_net = e_leaf * (ff_sky * e_sky * sigma * temp_sky ** 4 +
                                             e_leaf * longwave_grain_from_leaves +
                                             ff_soil * e_soil * sigma * temp_soil ** 4) \
                                   - 2 * e_leaf * sigma * t_leaf ** 4
                    latent_heat_loss = -lambda_ * evap
                    sensible_heat_net = -gb_h * (t_leaf - temp_air)
                    energy_balance = shortwave_abs + longwave_net + latent_heat_loss + sensible_heat_net
                    return energy_balance

                t_leaf0 = utils.kelvin_to_celsius(
                    optimize.newton_krylov(_VineEnergyX, utils.celsius_to_kelvin(t_leaf)))

                t_dict[vid] = t_leaf0

            t_new = t_dict

            # evaluation of leaf temperature conversion criterion
            error_dict = {vtx: abs(t_prev[vtx] - t_new[vtx]) for vtx in leaves}

            t_error = max(error_dict.values())
            t_error_trace.append(t_error)

            if t_error < t_error_crit:
                break
            else:
                try:
                    if abs(t_error_trace[-1] - t_error_trace[-2]) < t_error_crit:
                        it_step = max(0.01, it_step / 2.)
                except IndexError:
                    pass

                t_next = {}
                for vtx_id in t_new.keys():
                    tx = t_prev[vtx_id] + it_step * (t_new[vtx_id] - t_prev[vtx_id])
                    t_next[vtx_id] = tx

                t_prev = t_next

    # matrix iterative calculation of leaves temperature ('not solo' case)
    else:
        it = 1
        t_lst = []
        t_dict = {vid: Symbol('t%d' % vid) for vid in leaves}

        eq_lst = []
        t_leaf_lst = []
        for vid in leaves:
            shortwave_inc = properties['ei'][vid] / (0.48 * 4.6)  # Ei not Eabs
            ff_sky = properties['k_sky'][vid]
            ff_leaves = properties['k_leaves'][vid]
            ff_soil = properties['k_soil'][vid]
            gb_h = properties['gbh'][vid]
            evap = properties['ev'][vid]
            t_leaf = t_prev[vid]

            t_leaf_lst.append(t_leaf)
            t_lst.append(t_dict[vid])

            eq_aux = 0.
            for ivid in ff_leaves:
                if not g.node(ivid).label.startswith('soil'):
                    eq_aux += -ff_leaves[ivid] * ((t_dict[ivid]) ** 4)

            eq = (a_glob * shortwave_inc +
                  e_leaf * sigma * (ff_sky * e_sky * (temp_sky ** 4) +
                                    e_leaf * eq_aux + ff_soil * e_soil * (temp_sky ** 4) -
                                    2 * (t_dict[vid]) ** 4) -
                  lambda_ * evap - gb_h * Cp * (t_dict[vid] - temp_air))

            eq_lst.append(eq)

        tt = time.time()
        t_leaf0_lst = nsolve(eq_lst, t_lst, t_leaf_lst, verify=False) - 273.15
        print ("---%s seconds ---" % (time.time() - tt))

        t_new = {}
        for ivid, vid in enumerate(leaves):
            t_new[vid] = float(t_leaf0_lst[ivid])
            ivid += 1

    return t_new, it


def soil_temperature(g, meteo, temp_sky_eff, soil_label_prefix='other'):
    """Computes soil temperature

    Args:
        g: a multiscale tree graph object
        meteo (DataFrame): forcing meteorological variables
        t_sky_eff (float): [°C] effective sky temperature
        soil_label_prefix(str): prefix of soil nodes

    Returns:
        (double): [°C] soil temperature

    Notes:
        Heat loss into deeper soil layers is not considered.

    """

    relative_humidity, atm_pressure, temp_air = [float(meteo[x]) for x in ('hs', 'Pa', 'Tac')]

    temp_sky_eff = utils.celsius_to_kelvin(temp_sky_eff)

    soil_nodes = [g.node(vid) for vid in g.property('geometry') if g.node(vid).label.startswith(soil_label_prefix)][0]
    temp_leaves = utils.celsius_to_kelvin(mean(g.property('Tlc').values()))

    shortwave_inc = soil_nodes.Ei / (0.48 * 4.6)  # Ei not Eabs
    temp_soil = soil_nodes.Tsoil if 'Tsoil' in soil_nodes.properties() else temp_air

    def _SoilEnergyX(temp_soil):
        shortwave_abs = (1 - 0.25) * shortwave_inc  # 0.25 is rough estimation of albedo of a bare soil
        longwave_net = e_soil * sigma * (1. * e_sky * temp_sky_eff ** 4 +
                                         1. * e_leaf * temp_leaves ** 4 -
                                         temp_soil ** 4)
        latent_heat = -lambda_ * 0.06 * utils.vapor_pressure_deficit(temp_air, temp_soil,
                                                                     relative_humidity) / atm_pressure  # 0.06 is gM from Bailey 2016 AFM 218-219:146-160
        sensible_heat = -0.5 * Cp * utils.celsius_to_kelvin(
            temp_soil - temp_air)  # 0.5 is gH from Bailey 2016 AFM 218-219:146-160
        energy_balance = shortwave_abs + longwave_net + latent_heat + sensible_heat
        return energy_balance

    t_soil0 = utils.kelvin_to_celsius(
        optimize.newton_krylov(_SoilEnergyX, utils.celsius_to_kelvin(temp_soil)))

    soil_nodes.Tsoil = t_soil0

    return t_soil0


def forced_soil_temperature(meteo):
    """A very simple model of soil temperature

    Args:
        meteo (DataFrame): forcing meteorological variables

    Returns:
        (double): [°C] soil temperature

    """

    dt_soil = [3, 3, 3, 3, 3, 3, 3, 3, 10, 15, 20, 20, 20, 20, 20, 15, 6, 5, 4, 3, 3, 3, 3, 3]
    t_soil = meteo.Tac[0] + dt_soil[meteo.index.hour[0]]
    return t_soil