# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Hydraulic structure module of HydroShoot.

This module computes xylem water potential value at each node of the shoot.
"""

from copy import deepcopy
from math import exp, pi, log

import openalea.mtg.traversal as traversal
from numpy import array
from openalea.plantgl.all import surface as surf
from scipy import optimize

import hydroshoot.constants as cst
from hydroshoot.soil import calc_collar_water_potential


def conductivity_max(diameter, a=2.8, b=0.1, min_kmax=0.):
    """Computes the maximum conductivity of a hydraulic segment following Tyree and Zimermmann (2002)

    Args:
        diameter (float): [m] average diameter of the hydraulic segment
        a (float): [kg s-1 MPa-1] slope parameter between segment diameter and conductivity (ranges between 2.5 and
            2.8 according to Tyree and Zimmermann, 2002)
        b (float): [-] exponent to the segment diameter (ranges between 2.0 and 5.0 according to Tyree and
            Zimmermann, 2002)
        min_kmax (float): [kg s-1 m MPa-1] minimum value for the maximum conductivity

    Returns:
        (float): [kg s-1 m MPa-1] maximum conductivity of the hydraulic segment

    References:
        Tyree M, Zimmermann M. 2002.
            Xylem structure and the ascent of sap, Springer Series in Wood Science.
            (P. 147)
    """

    return max(min_kmax, float(a * (diameter ** b)))


def cavitation_factor(psi, model='tuzet', fifty_cent=-0.51, sig_slope=3):
    """Computes the effect of cavitation on xylem conductivity (i.e. the ratio of actual to maximum hydraulic
    conductivity).

    Args:
        psi (float): [MPa] water potential of the hydraulic segment
        model (str): one of 'misson' (logistic function with polynomial formula), 'tuzet' (logistic function with
            exponential formula), or 'linear' for linear reduction
        fifty_cent (float): [MPa] water potential at which the conductivity of the hydraulic segment drops to 50%
            of its maximum value
        sig_slope (float): a shape parameter controlling the slope of the S-curve (used only for 'misson' [-] or
            'tuzet' [MPa-1] models)

    Returns:
        (float): [-] the ratio of actual to maximum stem conductance (between 0 and 1)

    """

    if model == 'misson':
        k_reduction = 1. / (1. + (psi / fifty_cent) ** sig_slope)
    elif model == 'tuzet':
        k_reduction = (1. + exp(sig_slope * fifty_cent)) / (1. + exp(float(sig_slope * (fifty_cent - psi))))
    elif model == 'linear':
        k_reduction = 1 - min(0.95, float(psi / fifty_cent))
    else:
        raise ValueError("The 'model' argument must be one of the following ('misson','tuzet', 'linear').")

    return k_reduction


def def_param_soil(custom=None):
    """
    Returns a dictionary of classes of default soil hydrodynamic parameters for the model of van Genuchten-Muallem.
    For each soil class (`dict.key`), the data are organized as follows:
    name : (theta_r, theta_s, alpha[cm-1], n, k_sat[cm d-1])

    Args:
        custom (tuple): set of soil parameters ordered as mentioned above.

    References
        Carsel R., Parrish R., 1988.
            Developing joint probability distributions of soil water retention characteristics.
            Water Resources Research 24,755 – 769.
    """

    def_dict = {'Sand': (0.045, 0.430, 0.145, 2.68, 712.8),
                'Loamy_Sand': (0.057, 0.410, 0.124, 2.28, 350.2),
                'Sandy_Loam': (0.065, 0.410, 0.075, 1.89, 106.1),
                'Loam': (0.078, 0.430, 0.036, 1.56, 24.96),
                'Silt': (0.034, 0.460, 0.016, 1.37, 6.00),
                'Silty_Loam': (0.067, 0.450, 0.020, 1.41, 10.80),
                'Sandy_Clay_Loam': (0.100, 0.390, 0.059, 1.48, 31.44),
                'Clay_Loam': (0.095, 0.410, 0.019, 1.31, 6.24),
                'Silty_Clay_Loam': (0.089, 0.430, 0.010, 1.23, 1.68),
                'Sandy_Clay': (0.100, 0.380, 0.027, 1.23, 2.88),
                'Silty_Clay': (0.070, 0.360, 0.005, 1.09, 0.48),
                'Clay': (0.068, 0.380, 0.008, 1.09, 4.80)}

    if custom:
        def_dict['Custom'] = custom

    return def_dict


def k_soil_soil(psi, soil_class):
    """Gives the actual soil hydraulic conductivity following van Genuchten (1980)

    Args:
        psi (float): [MPa] bulk soil water potential
        soil_class (str): one of the soil classes proposed by Carsel and Parrish (1988), see :func:`def_param_soil` for
            details

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
    param = def_param_soil()[soil_class]
    theta_r, theta_s, alpha, n, k_sat = [param[i] for i in range(5)]
    m = 1. - 1. / n
    effective_saturation = 1. / (1. + abs(alpha * psi) ** n) ** m
    k_soil = k_sat * effective_saturation ** 0.5 * (1. - (1. - effective_saturation ** (1. / m)) ** m) ** 2

    return k_soil


def k_soil_root(k_soil, root_spacing, root_radius):
    """Computes the water conductivity at the soil-root interface according to Gardner (1960)

    Args:
        k_soil (float): [cm d-1] soil hydraulic conductivity
        root_spacing (float): [m] mean spacing between the neighbouring roots
        root_radius (float): [m] mean root radius

    Returns:
        (float): [cm d-1] water conductivity at the soil-root interface

    References:
        Gardner (1960)
            Dynamic aspects of water availability to plants.
            Soil science 89, 63 - 73.

    """
    return 4. * pi * k_soil / log((root_spacing / root_radius) ** 2)


def soil_water_potential(psi_soil_init, water_withdrawal, soil_class, soil_total_volume, psi_min=-3.):
    """Computes soil water potential following van Genuchten (1980)

    Args:
        psi_soil_init (float): [MPa] initial soil water potential
        water_withdrawal (float): [Kg T-1] water volume that is withdrawn from the soil (by transpiration for instance)
            during a time lapse T
        soil_class (str): one of the soil classes proposed by Carsel and Parrish (1988), see :func:`def_param_soil` for
            details
        soil_total_volume (float): [m3] total apparent volume of the soil (including solid, liquid and gaseous
            fractions)
        psi_min (float): [MPa] minimum allowable water potential

    Returns:
        (float): [MPa] soil water potential

    Notes:
        Strictly speaking, :arg:`psi_min` expresses rather the minimum water potential at the base of the plant shoot.

    References:
        van Genuchten M., 1980.
            A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
            Soil Science Society of America Journal 44, 892897.
    """

    psi_soil_init = min(-1e-6, psi_soil_init)
    theta_r, theta_s, alpha, n, k_sat = def_param_soil()[soil_class]

    m = 1. - 1. / n

    psi = psi_soil_init * 1.e6 / (cst.water_density * cst.gravitational_acceleration) * 100.  # MPa -> cm_H20
    theta_init = theta_r + (theta_s - theta_r) / (1. + abs(alpha * psi) ** n) ** m

    flux = water_withdrawal / cst.water_density  # kg T-1 -> m3 T-1

    porosity_volume = soil_total_volume * theta_s

    delta_theta = flux / porosity_volume  # [m3 m-3]

    theta = max(theta_r, theta_init - delta_theta)

    if theta <= theta_r:
        psi_soil = psi_min
    elif theta >= theta_s:
        psi_soil = 0.0
    else:
        def _water_retention(x):
            return 1. / (1. + abs(alpha * x) ** n) ** m - (theta - theta_r) / (theta_s - theta_r)

        psi_soil = optimize.fsolve(_water_retention, array(psi))[0] / (
                1.e6 / (cst.water_density * cst.gravitational_acceleration) * 100)

    return max(psi_min, float(psi_soil))


def hydraulic_prop(g, length_conv=1.e-2, a=2.6, b=2.0, min_kmax=0.):
    """Computes water flux `Flux` and maximum hydraulic conductivity `Kmax` of each hydraulic segment. Both properties
        are then attached to the corresponding mtg nodes.

    Args:
        g (openalea.mtg.MTG): a multiscale tree graph object
        length_conv (float): conversion coefficient from the length unit of the mtg to that of [1 m]
        a (float): [kg s-1 MPa-1] slope of the Kh(D) relationship, see :func:`conductivity_max` for details
        b (float): [-] exponent of the Kh(D) relationship, see :func:`conductivity_max` for details
        min_kmax (float): [kg s-1 m MPa-1] minimum value for the maximum conductivity, see :func:`conductivity_max`
            for details

    Returns:
        (openalea.mtg.MTG): the multiscale tree graph object

    Notes:
        The units of **Flux** and **Kmax** properties are related:
            - if `Flux` is given as water flux [kg s-1], then `Kmax` is in [kg m s-1 Pa-1] (or [kg m s-1 MPa-1])
            - else if `Flux` is given as water flux density [kg m-2 s-1], then `Kmax` is in [kg m-1 s-1 Pa-1]
            (or [kg m-1 s-1 MPa-1])
        Transpiration flux density per leaf surface area `E` (propery of the :arg:`mtg`) must be in [mol m-2 s-1],
            otherwise, the `cst.water_molar_mass` value must be re-adapted

        The resulting water potential, calculated in :func:`transient_xylem_water_potential` is then given in [MPa]

    """

    vid_base = g.node(g.root).vid_base

    for vtx_id in traversal.post_order2(g, vid_base):
        n = g.node(vtx_id)
        if n.label.startswith('LI'):
            try:
                leaf_area = n.leaf_area * 1.
            except (AttributeError, TypeError):
                leaf_area = surf(n.geometry) * length_conv ** 2  # [m2]
                # Note: The surface of the leaf mesh is overestimated compared to allometry results
                # leaf_area = (0.0175*(n.Length*10.)**1.9057)*LengthConv**2 #[m2]
                n.leaf_area = leaf_area

            n.Flux = (n.E * cst.water_molar_mass * 1.e-3) * leaf_area  # [kg(H2O) s-1]
            # n.FluxC = ((n.An)*44.0095*1.e-9)*leaf_area # [kgCO2 s-1]
            n.FluxC = n.An * leaf_area  # [umol s-1]

        elif n.label.startswith(('in', 'cx', 'Pet')):
            n.Flux = sum([vtx.Flux for vtx in n.children()])
            diam = 0.5 * (n.TopDiameter + n.BotDiameter) * length_conv
            n.Kmax = conductivity_max(diam, a, b, min_kmax)

            n.FluxC = sum([vtx.FluxC for vtx in n.children()])

        elif n.label.startswith('rhyzo'):
            n.Flux = sum([vtx.Flux for vtx in n.children()])
            n.FluxC = sum([vtx.FluxC for vtx in n.children()])

            n.Kmax = None

    return g


def transient_xylem_water_potential(g, model='tuzet', length_conv=1.e-2, psi_soil=-0.6, psi_min=-3., fifty_cent=-0.51,
                                    sig_slope=1., root_spacing=0.013, root_radius=.0001,
                                    negligible_shoot_resistance=False,
                                    start_vid=None, stop_vid=None):
    """Computes a transient hydraulic structure of a plant shoot based on constant values of the hydraulic segments'
    conductivities. The hydraulic segments are assumed isotropic having only axial conductivities.

    Args:
        g (openalea.mtg.MTG): a multiscale tree graph object
        model (str): one of 'misson', 'tuzet' or 'linear', see :func:`cavitation_factor` for details
        length_conv (float): conversion coefficient from the length unit of the mtg to that of [1 m]
        psi_soil (float): [MPa] soil water potential
        psi_min (float): [MPa] minimum allowable water potential in the hydraulic segments
        fifty_cent (float): [MPa] water potential at which the conductivity of the hydraulic segment drops to 50%
            of its maximum value, see :func:`cavitation_factor` for details
        sig_slope (float): a shape parameter controlling the slope of the S-curve (used only for 'misson' [-] or
            'tuzet' [MPa-1] models), see :func:`cavitation_factor` for details
        root_spacing (float): [m] mean spacing between the neighbouring roots
        root_radius (float): [m] mean root radius
        negligible_shoot_resistance (bool): to consider (True) or not to consider (False) shoot resistance to xylem
            flow
        start_vid (int): vertex id from which the iteration starts (if `None` it is then taken to the basal element)
        stop_vid (int): vertex id at which the iteration breaks (if `None` iteration will include all elements up until
            the leaves)

    Notes:
        The flux through the stem must be given in [kg s-1]
        The hydraulic conductivity of the stem must be given in [kg m s-1 MPa-1]
    """

    vid_base = g.node(g.root).vid_base

    if start_vid is None:
        start_vid = vid_base

    for vtx_id in traversal.pre_order2(g, start_vid):
        if vtx_id == stop_vid:
            break
        else:
            n = g.node(vtx_id)
            p = n.parent()

            if n.label.startswith('LI'):
                n.psi_head = p.psi_head
            else:
                flux = n.properties()['Flux']
                length = n.properties()['Length'] * length_conv
                z_head = n.properties()['TopPosition'][2] * length_conv
                z_base = n.properties()['BotPosition'][2] * length_conv

                if vtx_id == vid_base:
                    # psi_base = psi_soil
                    psi_base = calc_collar_water_potential(
                        transpiration=flux,
                        bulk_soil_water_potential=psi_soil,
                        root_depth=1.2,
                        soil_class='Loam',
                        root_radius=0.0001,
                        root_length_density=2000.)
                else:
                    psi_base = p.psi_head

                try:
                    psi_head = n.psi_head
                    assert n.psi_head, "Water potential has a `None` value"
                except AttributeError:
                    psi_head = psi_base

                psi = 0.5 * (psi_head + psi_base)

                if n.label.startswith('rhyzo'):
                    cyl_diameter = n.TopDiameter * length_conv
                    depth = n.depth * length_conv
                    flux *= 8640. / (pi * cyl_diameter * depth)  # [cm d-1]

                    if n.label.startswith('rhyzo0'):
                        k_soil = k_soil_soil(psi, n.soil_class)  # [cm d-1]
                        g_act = k_soil_root(k_soil, root_spacing, root_radius)  # [cm d-1 m-1]
                        psi_head = max(psi_min, psi_base - (
                                flux / g_act) * cst.water_density * cst.gravitational_acceleration * 1.e-6)
                        k_act = None
                    else:
                        k_act = k_soil_soil(psi, n.soil_class)  # [cm d-1]
                        psi_head = max(psi_min, psi_base - (
                                length * flux / k_act) * cst.water_density * cst.gravitational_acceleration * 1.e-6)

                else:
                    k_max = n.properties()['Kmax']

                    if not negligible_shoot_resistance:
                        k_act = k_max * cavitation_factor(psi, model, fifty_cent, sig_slope)
                        psi_head = max(psi_min,
                                       psi_base - length * flux / k_act - (
                                               cst.water_density * cst.gravitational_acceleration * (
                                               z_head - z_base)) * 1.e-6)
                    else:
                        k_act = None
                        psi_head = max(psi_min, psi_base - (
                                cst.water_density * cst.gravitational_acceleration * (z_head - z_base)) * 1.e-6)

                n.psi_head = psi_head
                n.KL = k_act

    return


def xylem_water_potential(g, psi_soil=-0.8, model='tuzet', psi_min=-3.0, psi_error_crit=0.001, max_iter=100,
                          length_conv=1.E-2, fifty_cent=-0.51, sig_slope=0.1, root_spacing=0.013, root_radius=.0001,
                          negligible_shoot_resistance=False, start_vid=None, stop_vid=None, psi_step=0.5):
    """Computes the hydraulic structure of plant's shoot.

    Args:
        g (openalea.mtg.MTG): a multiscale tree graph object
        psi_soil (float): [MPa] soil water potential
        model (str): one of 'misson', 'tuzet' or 'linear', see :func:`cavitation_factor` for details
        psi_min (float): [MPa] minimum allowable water potential in the hydraulic segments
        psi_error_crit (float): [MPa] water potential difference threshold below which iterations cease
        max_iter (int): maximum number of iterations
        length_conv (float): conversion coefficient from the length unit of the mtg to that of [1 m]
        fifty_cent (float): [MPa] water potential at which the conductivity of the hydraulic segment drops to 50%
            of its maximum value, see :func:`cavitation_factor` for details
        sig_slope (float): a shape parameter controlling the slope of the S-curve (used only for 'misson' [-] or
            'tuzet' [MPa-1] models), see :func:`cavitation_factor` for details
        root_spacing (float): [m] mean spacing between the neighbouring roots
        root_radius (float): [m] mean root radius
        negligible_shoot_resistance (bool): to consider (True) or not to consider (False) shoot resistance to xylem
            flow
        start_vid (int or None): vertex id from which the iteration starts (if `None` it is then taken to the basal
            element)
        stop_vid (int or None): vertex id at which the iteration breaks (if `None` iteration will include all elements
            up until the leaves)
        psi_step (float): [m] reduction factor to the xylem water potential step between two consecutive iterations
            (between 0 and 1)

    Returns:
        (int): the number of iterations

    """

    counter = 0

    psi_error = psi_error_crit
    psi_error_trace = []
    ipsi_step = psi_step

    while psi_error >= psi_error_crit:
        psi_error = 0.

        psi_prev = deepcopy(g.property('psi_head'))

        transient_xylem_water_potential(g, model, length_conv, psi_soil, psi_min, fifty_cent, sig_slope,
                                        root_spacing, root_radius, negligible_shoot_resistance, start_vid, stop_vid)

        psi_new = deepcopy(g.property('psi_head'))

        psi_error_trace.append(psi_error)

        if counter > max_iter:
            "The numerical solution of the hydraulic structure did not converge."
        else:
            try:
                if abs(psi_error_trace[-1] - psi_error_trace[-2]) < psi_error_crit:
                    ipsi_step = max(0.01, ipsi_step / 2.)
            except IndexError:
                pass

            psi_avg_dict = {}

            for vtx_id in g.property('psi_head').keys():
                psi_error += abs(psi_prev[vtx_id] - psi_new[vtx_id])
                psi_avg_dict[vtx_id] = psi_prev[vtx_id] + psi_step * (psi_new[vtx_id] - psi_prev[vtx_id])

            g.properties()['psi_head'] = psi_avg_dict

        counter += 1

    return counter
