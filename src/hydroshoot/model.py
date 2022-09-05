# -*- coding: utf-8 -*-
"""This module performs a complete comutation scheme: irradiance absorption, gas-exchange, hydraulic structure,
energy-exchange, and soil water depletion, for each given time step.
"""
from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import openalea.mtg.traversal as traversal
from openalea.plantgl.all import Scene, surface
from pandas import DataFrame, date_range, DatetimeIndex

from hydroshoot import (initialisation, architecture, irradiance, exchange, hydraulic, energy, display, solver, utilities)
from hydroshoot.params import Params


def run(g, wd, scene=None, write_result=True, **kwargs):
    """
    Calculates leaf gas and energy exchange in addition to the hydraulic structure of an individual plant.

    :Parameters:
    - **g**: a multiscale tree graph object
    - **wd**: string, working directory
    - **scene**: PlantGl scene
    - **kwargs** can include:
        - **psi_soil**: [MPa] predawn soil water potential
        - **gdd_since_budbreak**: [°Cd] growing degree-day since bubreak
        - **sun2scene**: PlantGl scene, when prodivided, a sun object (sphere) is added to it
        - **soil_size**: [cm] length of squared mesh size
    """
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('+ Project: ', wd)
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    time_on = datetime.now()

    # Read user parameters
    params = Params(f'{wd}params.json')

    # ==============================================================================
    # Initialisation
    # ==============================================================================
    #   Weather data
    meteo_tab = initialisation.init_weather(path_project=wd, params=params)
    meteo = meteo_tab.loc[params.simulation.date_range, :]

    #   Determination of the simulation period
    sdate = params.simulation.date_beg

    # Time conversion
    time_conv = params.simulation.time_conv

    # Unit length conversion (from scene unit to the standard [m]) unit)
    unit_scene_length = params.simulation.unit_scene_length
    length_conv = params.simulation.length_conv

    # Reading available pre-dawn soil water potential data
    psi_pd = initialisation.init_soil_predawn_water_potential(path_project=wd, params=params, **kwargs)

    # Determination of cumulative degree-days parameter
    t_base = params.phenology.t_base
    budbreak_date = datetime.strptime(params.phenology.emdate, "%Y-%m-%d %H:%M:%S")

    if 'gdd_since_budbreak' in kwargs:
        gdd_since_budbreak = kwargs['gdd_since_budbreak']
    elif min(meteo_tab.index) <= budbreak_date:
        tdays = date_range(budbreak_date, sdate, freq='D')
        tmeteo = meteo_tab.loc[tdays, :].Tac.to_frame()
        tmeteo = tmeteo.set_index(DatetimeIndex(tmeteo.index).normalize())
        df_min = tmeteo.groupby(tmeteo.index).aggregate(np.min).Tac
        df_max = tmeteo.groupby(tmeteo.index).aggregate(np.max).Tac
        # df_tt = merge(df_max, df_min, how='inner', left_index=True, right_index=True)
        # df_tt.columns = ('max', 'min')
        # df_tt['gdd'] = df_tt.apply(lambda x: 0.5 * (x['max'] + x['min']) - t_base)
        # gdd_since_budbreak = df_tt['gdd'].cumsum()[-1]
        df_tt = 0.5 * (df_min + df_max)
        df_tt = df_tt.apply(lambda x: utilities.calc_effective_daily_temperature(
            temperature_air=x, temperature_base=t_base))
        gdd_since_budbreak = df_tt.cumsum()[-1]
    else:
        raise ValueError('Cumulative degree-days temperature is not provided.')

    print('GDD since budbreak = %d °Cd' % gdd_since_budbreak)

    # Determination of perennial structure arms (for grapevine)
    # arm_vid = {g.node(vid).label: g.node(vid).components()[0]._vid for vid in g.VtxList(Scale=2) if
    #            g.node(vid).label.startswith('arm')}

    # Soil reservoir dimensions (inter row, intra row, depth) [m]
    soil_dimensions = params.soil.soil_dimensions
    soil_total_volume = soil_dimensions[0] * soil_dimensions[1] * soil_dimensions[2]
    rhyzo_coeff = params.soil.rhyzo_coeff
    rhyzo_total_volume = rhyzo_coeff * np.pi * min(soil_dimensions[:2]) ** 2 / 4. * soil_dimensions[2]

    # Counter clockwise angle between the default X-axis direction (South) and
    # the real direction of X-axis.
    scene_rotation = params.irradiance.scene_rotation

    # Sky and cloud temperature [degreeC]
    t_sky = params.energy.t_sky
    t_cloud = params.energy.t_cloud

    # Topological location
    latitude = params.simulation.latitude
    longitude = params.simulation.longitude
    elevation = params.simulation.elevation
    geo_location = (latitude, longitude, elevation)

    # Pattern
    ymax, xmax = map(lambda dim: dim / length_conv, soil_dimensions[:2])
    pattern = ((-xmax / 2.0, -ymax / 2.0), (xmax / 2.0, ymax / 2.0))

    # Label prefix of the collar internode
    vtx_label = params.mtg_api.collar_label

    # Label prefix of the leaves
    leaf_lbl_prefix = params.mtg_api.leaf_lbl_prefix

    # Label prefices of stem elements
    stem_lbl_prefix = params.mtg_api.stem_lbl_prefix

    E_type = params.irradiance.E_type
    tzone = params.simulation.tzone
    turtle_sectors = params.irradiance.turtle_sectors
    icosphere_level = params.irradiance.icosphere_level
    turtle_format = params.irradiance.turtle_format

    energy_budget = params.simulation.energy_budget
    print('Energy_budget: %s' % energy_budget)

    # Optical properties
    opt_prop = params.irradiance.opt_prop

    print('Hydraulic structure: %s' % params.simulation.hydraulic_structure)

    psi_min = params.hydraulic.psi_min

    # Parameters of leaf Nitrogen content-related models
    Na_dict = params.exchange.Na_dict

    # Computation of the form factor matrix
    if energy_budget:
        print('Computing form factors...')
        g = energy.set_form_factors_simplified(
            g, pattern=pattern, infinite=True, leaf_lbl_prefix=leaf_lbl_prefix, turtle_sectors=turtle_sectors,
            icosphere_level=icosphere_level, unit_scene_length=unit_scene_length)

    # Soil class
    soil_class = params.soil.soil_class
    print('Soil class: %s' % soil_class)

    # Add rhyzosphere elements to mtg
    rhyzo_solution = params.soil.rhyzo_solution
    print('rhyzo_solution: %s' % rhyzo_solution)

    if rhyzo_solution:
        if not any(item.startswith('rhyzo') for item in g.property('label').values()):
            vid_collar = architecture.mtg_base(g, vtx_label=vtx_label)
            vid_base = architecture.add_soil_components(g, params.soil.rhyzo_radii, soil_dimensions, soil_class, vtx_label, 1. / length_conv)
        else:
            vid_collar = g.node(g.root).vid_collar
            vid_base = g.node(g.root).vid_base

            radius_prev = 0.

            for ivid, vid in enumerate(g.Ancestors(vid_collar)[1:]):
                radius = params.soil.rhyzo_radii[ivid] / length_conv
                g.node(vid).Length = radius - radius_prev
                g.node(vid).depth = soil_dimensions[2] / length_conv
                g.node(vid).TopDiameter = radius * 2.
                g.node(vid).BotDiameter = radius * 2.
                g.node(vid).soil_class = soil_class
                radius_prev = radius

    else:
        # Identifying and attaching the base node of a single MTG
        vid_collar = architecture.mtg_base(g, vtx_label=vtx_label)
        vid_base = vid_collar

    g.node(g.root).vid_base = vid_base
    g.node(g.root).vid_collar = vid_collar

    # Initializing sapflow to 0
    for vtx_id in traversal.pre_order2(g, g.node(g.root).vid_base):
        g.node(vtx_id).Flux = 0.

    # Addition of a soil element
    if 'Soil' not in g.properties()['label'].values():
        if 'soil_size' in kwargs:
            if kwargs['soil_size'] > 0.:
                architecture.add_soil(g, kwargs['soil_size'])
        else:
            architecture.add_soil(g, 500.)

    # Suppression of undesired geometry for light and energy calculations
    geom_prop = g.properties()['geometry']
    vidkeys = []
    for vid in g.properties()['geometry']:
        n = g.node(vid)
        if not n.label.startswith(('L', 'other', 'soil')):
            vidkeys.append(vid)
    [geom_prop.pop(x) for x in vidkeys]
    g.properties()['geometry'] = geom_prop

    # Attaching optical properties to MTG elements
    g = irradiance.optical_prop(g, leaf_lbl_prefix=leaf_lbl_prefix,
                                stem_lbl_prefix=stem_lbl_prefix, wave_band='SW',
                                opt_prop=opt_prop)

    # Estimation of Nitroen surface-based content according to Prieto et al. (2012)
    # Estimation of intercepted irradiance over past 10 days:
    if not 'Na' in g.property_names():
        print('Computing Nitrogen profile...')
        assert (sdate - min(
            meteo_tab.index)).days >= 10, 'Meteorological data do not cover 10 days prior to simulation date.'

        ppfd10_date = sdate + timedelta(days=-10)
        ppfd10t = date_range(ppfd10_date, sdate, freq='H')
        ppfd10_meteo = meteo_tab.loc[ppfd10t, :]
        caribu_source, RdRsH_ratio = irradiance.irradiance_distribution(ppfd10_meteo, geo_location, E_type,
                                                                        tzone, turtle_sectors, turtle_format,
                                                                        None, scene_rotation, None)

        # Compute irradiance interception and absorbtion
        g, caribu_scene = irradiance.hsCaribu(mtg=g,
                                              unit_scene_length=unit_scene_length,
                                              source=caribu_source, direct=False,
                                              infinite=True, nz=50, ds=0.5,
                                              pattern=pattern)

        g.properties()['Ei10'] = {vid: g.node(vid).Ei * time_conv / 10. / 1.e6 for vid in g.property('Ei').keys()}

        # Estimation of leaf surface-based nitrogen content:
        for vid in g.VtxList(Scale=3):
            if g.node(vid).label.startswith(leaf_lbl_prefix):
                g.node(vid).Na = exchange.leaf_Na(gdd_since_budbreak, g.node(vid).Ei10,
                                                  Na_dict['aN'],
                                                  Na_dict['bN'],
                                                  Na_dict['aM'],
                                                  Na_dict['bM'])

    # Define path to folder
    output_path = f'{wd}output{params.simulation.output_index}/'

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
    psi_stem = {}
    Tlc_dict = {}
    Ei_dict = {}
    an_dict = {}
    gs_dict = {}

    # The time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for date in meteo.time:
        print("=" * 72)
        print('Date', date, '\n')

        # Select of meteo data
        imeteo = meteo[meteo.time == date]

        # initiate wind speed
        g.properties()['u'] = energy.set_wind_speed(g=g, meteo=imeteo, leaf_lbl_prefix=leaf_lbl_prefix)

        # Add a date index to g
        g.date = datetime.strftime(date, "%Y%m%d%H%M%S")

        # Read soil water potntial at midnight
        if 'psi_soil' in kwargs:
            psi_soil = kwargs['psi_soil']
        else:
            if date.hour == 0:
                try:
                    psi_soil_init = psi_pd.loc[date, :][0]
                    psi_soil = psi_soil_init
                except KeyError:
                    pass
            # Estimate soil water potntial evolution due to transpiration
            else:
                psi_soil = hydraulic.soil_water_potential(psi_soil,
                                                          g.node(vid_collar).Flux * time_conv,
                                                          soil_class, soil_total_volume, psi_min)

        if 'sun2scene' not in kwargs or not kwargs['sun2scene']:
            sun2scene = None
        elif kwargs['sun2scene']:
            sun2scene = display.visu(g, def_elmnt_color_dict=True, scene=Scene())

        # Compute irradiance distribution over the scene
        caribu_source, RdRsH_ratio = irradiance.irradiance_distribution(imeteo, geo_location, E_type, tzone,
                                                                        turtle_sectors, turtle_format, sun2scene,
                                                                        scene_rotation, None)

        # Compute irradiance interception and absorbtion
        g, caribu_scene = irradiance.hsCaribu(mtg=g,
                                              unit_scene_length=unit_scene_length,
                                              source=caribu_source, direct=False,
                                              infinite=True, nz=50, ds=0.5,
                                              pattern=pattern)

        # g.properties()['Ei'] = {vid: 1.2 * g.node(vid).Ei for vid in g.property('Ei').keys()}

        # Trace intercepted irradiance on each time step
        rg_ls.append(sum([g.node(vid).Ei / (0.48 * 4.6) * surface(g.node(vid).geometry) * (length_conv ** 2) \
                          for vid in g.property('geometry') if g.node(vid).label.startswith('L')]))

        # Hack forcing of soil temperture (model of soil temperature under development)
        t_soil = energy.forced_soil_temperature(imeteo)

        # Climatic data for energy balance module
        # TODO: Change the t_sky_eff formula (cf. Gliah et al., 2011, Heat and Mass Transfer, DOI: 10.1007/s00231-011-0780-1)
        t_sky_eff = RdRsH_ratio * t_cloud + (1 - RdRsH_ratio) * t_sky

        solver.solve_interactions(g=g, meteo=imeteo, psi_soil=psi_soil, t_soil=t_soil, t_sky_eff=t_sky_eff,
                                  length_conv=length_conv, time_conv=time_conv,
                                  rhyzo_total_volume=rhyzo_total_volume, params=params)

        # Write mtg to an external file
        if scene is not None:
            architecture.mtg_save(g, scene, output_path)

        # Plot stuff..
        sapflow.append(g.node(vid_collar).Flux)
        # sapEast.append(g.node(arm_vid['arm1']).Flux)
        # sapWest.append(g.node(arm_vid['arm2']).Flux)

        an_ls.append(g.node(vid_collar).FluxC)

        psi_stem[date] = deepcopy(g.property('psi_head'))
        Tlc_dict[date] = deepcopy(g.property('Tlc'))
        Ei_dict[date] = deepcopy(g.property('Eabs'))
        an_dict[date] = deepcopy(g.property('An'))
        gs_dict[date] = deepcopy(g.property('gs'))

        print('---------------------------')
        print('psi_soil', round(psi_soil, 4))
        print('psi_collar', round(g.node(3).psi_head, 4))
        print('psi_leaf', round(np.median([g.node(vid).psi_head for vid in g.property('gs').keys()]), 4))
        print('')
        # print('Rdiff/Rglob ', RdRsH_ratio)
        # print('t_sky_eff ', t_sky_eff)
        print('gs', np.median(list(g.property('gs').values())))
        print('flux H2O', round(g.node(vid_collar).Flux * 1000. * time_conv, 4))
        print('flux C2O', round(g.node(vid_collar).FluxC, 4))
        print('Tleaf ', round(np.median([g.node(vid).Tlc for vid in g.property('gs').keys()]), 2),
              'Tair ', round(imeteo.Tac[0], 4))
        print('')
        print("=" * 72)

    # End time loop +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Write output
    # Plant total transpiration
    sapflow = [flow * time_conv * 1000. for flow in sapflow]

    # sapEast, sapWest = [np.array(flow) * time_conv * 1000. for i, flow in enumerate((sapEast, sapWest))]

    # Median leaf temperature
    t_ls = [np.median(list(Tlc_dict[date].values())) for date in meteo.time]

    # Intercepted global radiation
    rg_ls = np.array(rg_ls) / (soil_dimensions[0] * soil_dimensions[1])

    results_dict = {
        'Rg': rg_ls,
        'An': an_ls,
        'E': sapflow,
        # 'sapEast': sapEast,
        # 'sapWest': sapWest,
        'Tleaf': t_ls
    }

    # Results DataFrame
    results_df = DataFrame(results_dict, index=meteo.time)

    # Write
    if write_result:
        results_df.to_csv(output_path + 'time_series.output',
                          sep=';', decimal='.')

    time_off = datetime.now()

    print("")
    print("beg time", time_on)
    print("end time", time_off)
    print("--- Total runtime: %d minute(s) ---" % int((time_off - time_on).seconds / 60.))

    return results_df
