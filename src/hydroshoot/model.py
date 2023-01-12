# -*- coding: utf-8 -*-
"""This module performs a complete comutation scheme: irradiance absorption, gas-exchange, hydraulic structure,
energy-exchange, and soil water depletion, for each given time step.
"""
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import numpy as np
from openalea.mtg.mtg import MTG
from openalea.plantgl.all import Scene, surface
from pandas import DataFrame

from hydroshoot import (architecture, solver, io)
from hydroshoot.energy import calc_effective_sky_temperature
from hydroshoot.initialisation import init_model, init_hourly


def run(g: MTG, wd: Path, scene: Scene = None, write_result: bool = True, path_output: Path = None,
        **kwargs) -> DataFrame:
    """Calculates leaf gas and energy exchange in addition to the hydraulic structure of an individual plant.

    Args:
        g: mtg object
        wd: working directory
        scene: PlantGl scene (default None)
        write_result: if True then hourly plant-scale outputs are written into a CSV file
        path_output: summary data output file path
        kwargs: can include:
            psi_soil (float): [MPa] predawn soil water potential
            gdd_since_budbreak (float): [°Cd] growing degree-day since bubreak
            sun2scene (Scene): PlantGl scene, when prodivided, a sun object (sphere) is added to it
            soil_size (float): [cm] length of squared mesh size
            leaf_nitrogen (dict): leaf nitrogen content per area (key=(int) mtg leaf vertex, value=(float) nitrogen content)
            leaf_ppfd (dict of dict): incident and absorbed PPFD by each leaf per each simulated hour
                key:(datetime) simulated datetime, value:
                    key:'Ei', value: (key: (int) mtg leaf vertex, value: (incident PPFD)),
                    key:'Eabs', value: (key: (int) mtg leaf vertex, value: (absorbed PPFD))
            form_factors (dict of dict): form factors for the sky, leaves and the soil
                key=(str) one of ('ff_sky', 'ff_leaves', 'ff_soil'), value=(key=(int) mtg leaf vertex, value=(form factor)

    Returns:
        Absorbed whole plant global irradiance (Rg), net photosynthesis (An), transpiration (E) and
            median leaf temperature (Tleaf)

    """
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('+ Project: ', wd)
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    time_on = datetime.now()

    # Read user parameters
    inputs = io.HydroShootInputs(
        path_project=wd,
        scene=scene,
        write_result=write_result,
        path_output_file=path_output,
        **kwargs)
    io.verify_inputs(g=g, inputs=inputs)
    params = inputs.params

    # ==============================================================================
    # Initialisation
    # ==============================================================================
    # Conversions
    time_conv = params.simulation.conv_to_second

    io.print_sim_infos(inputs=inputs)

    # Determination of perennial structure arms (for grapevine)
    # arm_vid = {g.node(vid).label: g.node(vid).components()[0]._vid for vid in g.VtxList(Scale=2) if
    #            g.node(vid).label.startswith('arm')}

    g = init_model(g=g, inputs=inputs)

    # Save geometry in an external file
    # HSArc.mtg_save_geometry(scene, output_path)

    # ==============================================================================
    # Simulations
    # ==============================================================================

    sapflow = []
    # sapEast = []
    # sapWest = []
    an_ls = []
    rg_ls = []
    leaf_temperature_dict = {}

    # The time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    inputs_hourly = io.HydroShootHourlyInputs(psi_soil=inputs.psi_soil_forced, sun2scene=inputs.sun2scene)

    for date in params.simulation.date_range:
        print("=" * 72)
        print(f'Date: {date}\n')

        # Select meteo data
        inputs_hourly.update(g=g, date_sim=date, hourly_weather=inputs.weather[inputs.weather.index == date],
                             psi_pd=inputs.psi_pd, params=params)

        g, diffuse_to_total_irradiance_ratio = init_hourly(
            g=g, inputs_hourly=inputs_hourly, leaf_ppfd=inputs.leaf_ppfd, params=params)

        inputs_hourly.sky_temperature = calc_effective_sky_temperature(
            diffuse_to_total_irradiance_ratio=diffuse_to_total_irradiance_ratio,
            temperature_cloud=params.energy.t_cloud,
            temperature_sky=params.energy.t_sky)

        solver.solve_interactions(
            g=g, meteo=inputs_hourly.weather.loc[date], psi_soil=inputs_hourly.psi_soil,
            t_soil=inputs_hourly.soil_temperature, t_sky_eff=inputs_hourly.sky_temperature, params=params)

        # Write mtg to an external file
        if scene is not None:
            architecture.save_mtg(g=g, scene=scene, file_path=inputs.path_output_dir)

        # Plot stuff..
        sapflow.append(g.node(g.node(g.root).vid_collar).Flux)
        # sapEast.append(g.node(arm_vid['arm1']).Flux)
        # sapWest.append(g.node(arm_vid['arm2']).Flux)

        # Trace intercepted irradiance on each time step
        rg_ls.append(
            sum([g.node(vid).Ei / (0.48 * 4.6) * surface(g.node(vid).geometry) * (params.simulation.conv_to_meter ** 2)
                 for vid in g.property('geometry') if g.node(vid).label.startswith('L')]))

        an_ls.append(g.node(g.node(g.root).vid_collar).FluxC)

        leaf_temperature_dict[date] = deepcopy(g.property('Tlc'))

        print('---------------------------')
        print(f'psi_soil {inputs_hourly.psi_soil:.4f}')
        print(f'psi_collar {g.node(g.node(g.root).vid_collar).psi_head:.4f}')
        print(f'psi_leaf {np.median([g.node(vid).psi_head for vid in g.property("gs").keys()]):.4f}')
        print('')
        # print('Rdiff/Rglob ', RdRsH_ratio)
        # print('t_sky_eff ', t_sky_eff)
        print(f'gs: {np.median(list(g.property("gs").values())):.4f}')
        print(f'flux H2O {g.node(g.node(g.root).vid_collar).Flux * 1000. * time_conv:.4f}')
        print(f'flux C2O {g.node(g.node(g.root).vid_collar).FluxC}')
        print(f'Tleaf {np.median([g.node(vid).Tlc for vid in g.property("gs").keys()]):.2f}', ' ',
              f'Tair {inputs_hourly.weather.loc[date, "Tac"]:.2f}')
        print('')
        print("=" * 72)

    # End time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Write output
    # Plant total transpiration
    sapflow = [flow * time_conv * 1000. for flow in sapflow]

    # sapEast, sapWest = [np.array(flow) * time_conv * 1000. for i, flow in enumerate((sapEast, sapWest))]

    # Median leaf temperature
    t_ls = [np.median(list(leaf_temperature_dict[date].values())) for date in params.simulation.date_range]

    # Intercepted global radiation
    rg_ls = np.array(rg_ls) / (params.planting.spacing_on_row * params.planting.spacing_between_rows)

    # Results DataFrame
    results_df = DataFrame({
        'Rg': rg_ls,
        'An': an_ls,
        'E': sapflow,
        # 'sapEast': sapEast,
        # 'sapWest': sapWest,
        'Tleaf': t_ls}, index=params.simulation.date_range)

    # Write
    if write_result:
        results_df.to_csv(inputs.path_output_file, sep=';', decimal='.')

    time_off = datetime.now()

    print("")
    print("beg time", time_on)
    print("end time", time_off)
    print(f"--- Total runtime: {(time_off - time_on).seconds} sec ---")
    return results_df
