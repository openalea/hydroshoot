# -*- coding: utf-8 -*-
"""
@author: Rami ALBASHA

Hydraulic structure module of HydroShoot.

This module computes xylem water potential value at each node of the shoot.
"""
from copy import deepcopy
from math import exp
from typing import Callable

import openalea.mtg.traversal as traversal
from openalea.plantgl.all import surface as surf

import openalea.hydroshoot.constants as cst


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


def transient_xylem_water_potential(g, calc_collar_water_potential: Callable,
                                    model='tuzet', length_conv=1.e-2, psi_soil=-0.6, psi_min=-3., fifty_cent=-0.51,
                                    sig_slope=1., negligible_shoot_resistance=False, start_vid=None, stop_vid=None):
    """Computes a transient hydraulic structure of a plant shoot based on constant values of the hydraulic segments'
    conductivities. The hydraulic segments are assumed isotropic having only axial conductivities.

    Args:
        g (openalea.mtg.MTG): a multiscale tree graph object
        calc_collar_water_potential: a function that takes 'transpiration' and 'soil_water_potential' for inputs and
            returns the collar water potential as scalar output
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
                    psi_base = calc_collar_water_potential(
                        transpiration=flux,
                        soil_water_potential=psi_soil)
                else:
                    psi_base = p.psi_head

                try:
                    psi_head = n.psi_head
                    assert n.psi_head is not None, "Water potential has a `None` value"
                except AttributeError:
                    psi_head = psi_base

                psi = 0.5 * (psi_head + psi_base)

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


def xylem_water_potential(g, calc_collar_water_potential,
                          psi_soil=-0.8, model='tuzet', psi_min=-3.0, psi_error_crit=0.001, max_iter=100,
                          length_conv=1.E-2, fifty_cent=-0.51, sig_slope=0.1, negligible_shoot_resistance=False,
                          start_vid=None, stop_vid=None, psi_step=0.5):
    """Computes the hydraulic structure of plant's shoot.

    Args:
        g (openalea.mtg.MTG): a multiscale tree graph object
        calc_collar_water_potential: a function that takes 'transpiration' and 'soil_water_potential' for inputs and
            returns the collar water potential as scalar output
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

        transient_xylem_water_potential(
            g=g,
            calc_collar_water_potential=calc_collar_water_potential,
            model=model,
            length_conv=length_conv,
            psi_soil=psi_soil,
            psi_min=psi_min,
            fifty_cent=fifty_cent,
            sig_slope=sig_slope,
            negligible_shoot_resistance=negligible_shoot_resistance,
            start_vid=start_vid,
            stop_vid=stop_vid)

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
