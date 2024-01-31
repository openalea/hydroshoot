from copy import deepcopy

import openalea.mtg.traversal as traversal

from openalea.hydroshoot import hydraulic, exchange, energy
from openalea.hydroshoot.architecture import get_leaves
from openalea.hydroshoot.soil import update_soil_water_potential


def solve_interactions(g, meteo, psi_soil, t_soil, t_sky_eff, params, calc_collar_water_potential):
    """Computes gas-exchange, energy and hydraulic structure of plant's shoot jointly.

    Args:
        g: MTG object
        meteo (DataFrame): forcing meteorological variables
        psi_soil (float): [MPa] soil (root zone) water potential
        t_soil (float): [degreeC] soil surface temperature
        t_sky_eff (float): [degreeC] effective sky temperature
        params (params): [-] :class:`hydroshoot.params.Params()` object
        calc_collar_water_potential: a function that takes 'transpiration' and 'soil_water_potential' for inputs and
            returns the collar water potential as scalar output

    """
    length_conv = params.simulation.conv_to_meter
    time_conv = params.simulation.conv_to_second

    is_energy_budget = params.simulation.is_energy_budget

    photosynthesis_params = params.exchange.par_photo
    stomatal_conductance_params = params.exchange.par_gs
    turbulence_resistance = params.exchange.rbt

    xylem_k_max = params.hydraulic.Kx_dict
    xylem_k_cavitation = params.hydraulic.par_K_vul
    psi_min = params.hydraulic.psi_min

    irradiance_type2 = params.irradiance.E_type2

    temp_step = params.numerical_resolution.t_step
    psi_step = params.numerical_resolution.psi_step
    max_iter = params.numerical_resolution.max_iter
    psi_error_threshold = params.numerical_resolution.psi_error_threshold
    temp_error_threshold = params.numerical_resolution.t_error_threshold

    modelx, psi_critx, slopex = [xylem_k_cavitation[ikey] for ikey in ('model', 'fifty_cent', 'sig_slope')]

    leaf_ids = get_leaves(g=g, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    # Initializations ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Initialize all xylem potential values to soil water potential
    for vtx_id in traversal.pre_order2(g, g.node(g.root).vid_base):
        g.node(vtx_id).psi_head = psi_soil

    # If leaf temperature to be calculated, calculate the boundary layer conductance to heat transfer
    if is_energy_budget:
        g = energy.calc_heat_boundary_layer_conductance(
            g=g, leaf_ids=leaf_ids, unit_scene_length=params.simulation.unit_scene_length)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Temperature loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    t_error_trace = []
    it_step = temp_step

    t_error = temp_error_threshold
    it = 0
    while (t_error >= temp_error_threshold) and (it < max_iter):
        t_prev = deepcopy(g.property('Tlc'))

        # Hydraulic loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if params.simulation.is_hydraulic_structure:
            psi_error_trace = []
            ipsi_step = psi_step
            for ipsi in range(max_iter):
                psi_prev = deepcopy(g.property('psi_head'))

                # Compute gas-exchange fluxes. Leaf T and Psi are from prev calc loop
                exchange.set_gas_exchange_rates(
                    g=g, photo_params=photosynthesis_params, gs_params=stomatal_conductance_params,
                    relative_humidity=meteo['hs'], air_co2=meteo['Ca'],
                    atmospheric_pressure=meteo['Pa'], E_type2=irradiance_type2, leaf_ids=leaf_ids,
                    rbt=turbulence_resistance)

                # Compute sap flow and hydraulic properties
                hydraulic.hydraulic_prop(g, length_conv=length_conv,
                                         a=xylem_k_max['a'], b=xylem_k_max['b'], min_kmax=xylem_k_max['min_kmax'])

                # Update soil water status
                psi_base = update_soil_water_potential(
                    psi_soil_init=psi_soil, water_withdrawal=g.node(g.node(g.root).vid_collar).Flux * time_conv,
                    soil_class=params.soil.soil_class, soil_total_volume=params.soil.soil_volume,
                    psi_min=psi_min)

                # Compute xylem water potential
                n_iter_psi = hydraulic.xylem_water_potential(
                    g=g, calc_collar_water_potential=calc_collar_water_potential,
                    psi_soil=psi_base, model=modelx, psi_min=psi_min, psi_error_crit=psi_error_threshold,
                    max_iter=max_iter, length_conv=length_conv, fifty_cent=psi_critx, sig_slope=slopex,
                    negligible_shoot_resistance=params.simulation.is_negligible_shoot_resistance,
                    start_vid=g.node(g.root).vid_base, stop_vid=None, psi_step=psi_step)

                psi_new = g.property('psi_head')

                # Evaluate xylem conversion criterion
                psi_error_dict = {}
                for vtx_id in g.property('psi_head').keys():
                    psi_error_dict[vtx_id] = abs(psi_prev[vtx_id] - psi_new[vtx_id])

                psi_error = max(psi_error_dict.values())
                psi_error_trace.append(psi_error)

                print('psi_error = ', round(psi_error, 3), ':: Nb_iter = %d' % n_iter_psi, 'ipsi_step = %f' % ipsi_step)

                # Manage xylem water potential step to ensure convergence
                if psi_error < psi_error_threshold:
                    break
                else:
                    try:
                        if psi_error_trace[-1] >= psi_error_trace[-2] - psi_error_threshold:
                            ipsi_step = max(0.05, ipsi_step / 2.)
                    except IndexError:
                        pass

                    psi_new_dict = {}
                    for vtx_id in psi_new.keys():
                        psix = psi_prev[vtx_id] + ipsi_step * (psi_new[vtx_id] - psi_prev[vtx_id])
                        psi_new_dict[vtx_id] = psix

                    g.properties()['psi_head'] = psi_new_dict

        else:
            # Compute gas-exchange fluxes. Leaf T and Psi are from prev calc loop
            exchange.set_gas_exchange_rates(
                g=g, photo_params=photosynthesis_params, gs_params=stomatal_conductance_params,
                relative_humidity=meteo['hs'], air_co2=meteo['Ca'],
                atmospheric_pressure=meteo['Pa'], E_type2=irradiance_type2, leaf_ids=leaf_ids,
                rbt=turbulence_resistance)

            # Compute sap flow and hydraulic properties
            hydraulic.hydraulic_prop(g, length_conv=length_conv,
                                     a=xylem_k_max['a'], b=xylem_k_max['b'], min_kmax=xylem_k_max['min_kmax'])

        # End Hydraulic loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Compute leaf temperature
        if is_energy_budget:
            it += 1
            if params.energy.solo:
                g.properties()['Tlc'], t_iter = energy.calc_leaf_temperature(
                    g=g, t_soil=t_soil, t_sky_eff=t_sky_eff, leaf_ids=leaf_ids,
                    max_iter=max_iter, t_error_crit=temp_error_threshold, t_step=temp_step)
            else:
                g.properties()['Tlc'], t_iter = energy.calc_leaf_temperature2(
                    g=g, t_soil=t_soil, t_sky_eff=t_sky_eff, leaf_ids=leaf_ids)

            # t_iter_list.append(t_iter)
            t_new = deepcopy(g.property('Tlc'))

            # Evaluation of leaf temperature conversion creterion
            error_dict = {vtx: abs(t_prev[vtx] - t_new[vtx]) for vtx in g.property('Tlc').keys()}

            t_error = round(max(error_dict.values()), 3)
            print('t_error = ', t_error, 'counter =', it, 't_iter = ', t_iter, 'it_step = ', it_step)
            t_error_trace.append(t_error)

            try:
                if t_error_trace[-1] >= t_error_trace[-2] - temp_error_threshold:
                    it_step = max(0.001, it_step / 2.)
            except IndexError:
                pass

            t_new_dict = {}
            for vtx_id in t_new.keys():
                tx = t_prev[vtx_id] + it_step * (t_new[vtx_id] - t_prev[vtx_id])
                t_new_dict[vtx_id] = tx

            g.properties()['Tlc'] = t_new_dict
            g = energy.set_local_vpd(g=g, relative_humidity=meteo['hs'], leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)
        else:
            t_error = 0

    # End temperature loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
