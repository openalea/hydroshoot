from copy import deepcopy

import openalea.mtg.traversal as traversal

from hydroshoot import hydraulic, exchange, energy


def solve_interactions(g, meteo, psi_soil, t_soil, t_sky_eff, params):
    """Computes gas-exchange, energy and hydraulic structure of plant's shoot jointly.

    Args:
        g: MTG object
        meteo (DataFrame): forcing meteorological variables
        psi_soil (float): [MPa] soil (root zone) water potential
        t_soil (float): [degreeC] soil surface temperature
        t_sky_eff (float): [degreeC] effective sky temperature
        params (params): [-] :class:`hydroshoot.params.Params()` object

    """
    length_conv = params.simulation.conv_to_meter
    time_conv = params.simulation.conv_to_second
    rhyzo_total_volume = params.soil.rhyzo_total_volume

    hydraulic_structure = params.simulation.hydraulic_structure
    negligible_shoot_resistance = params.simulation.negligible_shoot_resistance
    soil_water_deficit = params.simulation.soil_water_deficit
    energy_budget = params.simulation.energy_budget

    par_photo = params.exchange.par_photo
    par_gs = params.exchange.par_gs
    rbt = params.exchange.rbt

    xylem_k_max = params.hydraulic.Kx_dict
    xylem_k_cavitation = params.hydraulic.par_K_vul
    psi_min = params.hydraulic.psi_min

    solo = params.energy.solo

    irradiance_type2 = params.irradiance.E_type2

    leaf_lbl_prefix = params.mtg_api.leaf_lbl_prefix

    soil_class = params.soil.soil_class

    temp_step = params.numerical_resolution.t_step
    psi_step = params.numerical_resolution.psi_step
    max_iter = params.numerical_resolution.max_iter
    psi_error_threshold = params.numerical_resolution.psi_error_threshold
    temp_error_threshold = params.numerical_resolution.t_error_threshold

    modelx, psi_critx, slopex = [xylem_k_cavitation[ikey] for ikey in ('model', 'fifty_cent', 'sig_slope')]

    if hydraulic_structure:
        assert (par_gs['model'] != 'vpd'), "Stomatal conductance model should be linked to the hydraulic strucutre"
    else:
        par_gs['model'] = 'vpd'
        negligible_shoot_resistance = True

        print("par_gs: 'model' is forced to 'vpd'")
        print("negligible_shoot_resistance is forced to True.")

    # Initializations ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Initialize all xylem potential values to soil water potential
    for vtx_id in traversal.pre_order2(g, g.node(g.root).vid_base):
        g.node(vtx_id).psi_head = psi_soil

    # If leaf temperature to be calculated, calculate the boundary layer conductance to heat transfer
    if energy_budget:
        g = energy.calc_heat_boundary_layer_conductance(
            g=g, leaf_label_prefix=leaf_lbl_prefix, unit_scene_length=params.simulation.unit_scene_length)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Temperature loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    t_error_trace = []
    it_step = temp_step

    # Initialize leaf  temperature to air temperature
    g.properties()['Tlc'] = energy.set_leaf_temperature_to_air_temperature(g, meteo, leaf_lbl_prefix)

    for it in range(max_iter):
        t_prev = deepcopy(g.property('Tlc'))

        # Hydraulic loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if hydraulic_structure:
            psi_error_trace = []
            ipsi_step = psi_step
            for ipsi in range(max_iter):
                psi_prev = deepcopy(g.property('psi_head'))

                # Compute gas-exchange fluxes. Leaf T and Psi are from prev calc loop
                exchange.gas_exchange_rates(
                    g=g, photo_params=par_photo, gs_params=par_gs, air_temperature=meteo['Tac'],
                    relative_humidity=meteo['hs'], air_co2=meteo['Ca'], atmospheric_pressure=meteo['Pa'],
                    E_type2=irradiance_type2, leaf_lbl_prefix=leaf_lbl_prefix, rbt=rbt)

                # Compute sap flow and hydraulic properties
                hydraulic.hydraulic_prop(g, length_conv=length_conv,
                                         a=xylem_k_max['a'], b=xylem_k_max['b'], min_kmax=xylem_k_max['min_kmax'])

                # Update soil water status
                psi_base = hydraulic.soil_water_potential(
                    psi_soil_init=psi_soil, water_withdrawal=g.node(g.node(g.root).vid_collar).Flux * time_conv,
                    soil_class=soil_class, soil_total_volume=rhyzo_total_volume, psi_min=psi_min)

                if soil_water_deficit:
                    psi_base = max(-1.3, psi_base)
                else:
                    psi_base = max(-0.7, psi_base)

                # Compute xylem water potential
                n_iter_psi = hydraulic.xylem_water_potential(g, psi_soil=psi_base, model=modelx, psi_min=psi_min,
                                                             psi_error_crit=psi_error_threshold, max_iter=max_iter,
                                                             length_conv=length_conv, fifty_cent=psi_critx,
                                                             sig_slope=slopex,
                                                             root_spacing=params.soil.avg_root_spacing,
                                                             root_radius=params.soil.avg_root_radius,
                                                             negligible_shoot_resistance=negligible_shoot_resistance,
                                                             start_vid=g.node(g.root).vid_base, stop_vid=None,
                                                             psi_step=psi_step)

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
            exchange.gas_exchange_rates(
                g=g, photo_params=par_photo, gs_params=par_gs, air_temperature=meteo['Tac'],
                relative_humidity=meteo['hs'], air_co2=meteo['Ca'], atmospheric_pressure=meteo['Pa'],
                E_type2=irradiance_type2, leaf_lbl_prefix=leaf_lbl_prefix, rbt=rbt)

            # Compute sap flow and hydraulic properties
            hydraulic.hydraulic_prop(g, length_conv=length_conv,
                                     a=xylem_k_max['a'], b=xylem_k_max['b'], min_kmax=xylem_k_max['min_kmax'])

        # End Hydraulic loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Compute leaf temperature
        if energy_budget:
            g.properties()['Tlc'], t_iter = energy.leaf_temperature(g, meteo, t_soil, t_sky_eff,
                                                                    solo=solo,
                                                                    leaf_lbl_prefix=leaf_lbl_prefix, max_iter=max_iter,
                                                                    t_error_crit=temp_error_threshold, t_step=temp_step)

            # t_iter_list.append(t_iter)
            t_new = deepcopy(g.property('Tlc'))

            # Evaluation of leaf temperature conversion creterion
            error_dict = {vtx: abs(t_prev[vtx] - t_new[vtx]) for vtx in g.property('Tlc').keys()}

            t_error = round(max(error_dict.values()), 3)
            print('t_error = ', t_error, 'counter =', it, 't_iter = ', t_iter, 'it_step = ', it_step)
            t_error_trace.append(t_error)

            # Manage temperature step to ensure convergence
            if t_error < temp_error_threshold:
                break
            else:
                assert (it <= max_iter), 'The energy budget solution did not converge.'

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

    # End temperature loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
