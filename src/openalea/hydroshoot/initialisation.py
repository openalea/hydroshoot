from datetime import datetime, timedelta
from typing import Callable

import openalea.mtg.traversal as traversal
from openalea.mtg.mtg import MTG
from pandas import DataFrame
from pandas import date_range

from openalea.hydroshoot import soil
from openalea.hydroshoot.architecture import get_mtg_base, add_soil_surface_mesh, get_leaves
from openalea.hydroshoot.energy import (set_form_factors_simplified, set_local_wind_speed, set_local_air_temperature,
                               set_local_vpd, set_leaf_temperature_to_air_temperature)
from openalea.hydroshoot.exchange import leaf_Na
from openalea.hydroshoot.io import HydroShootInputs, HydroShootHourlyInputs
from openalea.hydroshoot.irradiance import irradiance_distribution, hsCaribu, set_optical_properties
from openalea.hydroshoot.params import Params
from openalea.hydroshoot.preprocess import calc_gdd_since_budbreak


def calc_nitrogen_distribution(g: MTG, gdd_since_budbreak: float, weather: DataFrame, params: Params) -> float:
    if gdd_since_budbreak is None:
        gdd_since_budbreak = calc_gdd_since_budbreak(
            weather=weather,
            date_beg_sim=params.simulation.date_beg,
            date_budbreak=params.phenology.date_budbreak,
            temperature_base=params.phenology.t_base)

    date_beg_sim = params.simulation.date_beg
    date_range_last_10_days = date_range(date_beg_sim + timedelta(days=-10), date_beg_sim, freq='H')
    caribu_source, _ = irradiance_distribution(
        meteo=weather.loc[date_range_last_10_days, :],
        geo_location=params.simulation.geo_location,
        irradiance_unit=params.irradiance.E_type,
        time_zone=params.simulation.tzone,
        turtle_sectors=params.irradiance.turtle_sectors,
        turtle_format=params.irradiance.turtle_format,
        sun2scene=None,
        rotation_angle=params.planting.scene_rotation)

    # Compute irradiance interception and absorption
    g, _ = hsCaribu(
        mtg=g,
        unit_scene_length=params.simulation.unit_scene_length,
        source=caribu_source,
        direct=False,
        infinite=True,
        nz=50,
        ds=0.5,
        pattern=params.irradiance.pattern)

    g.properties()['Ei10'] = {vid: g.node(vid).Ei * params.simulation.conv_to_second / 10. / 1.e6
                              for vid in g.property('Ei').keys()}
    g.properties().pop('Ei')

    # Estimation of leaf surface-based nitrogen content:
    for vid in get_leaves(g=g, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix):
        g.node(vid).Na = leaf_Na(
            age_gdd=gdd_since_budbreak,
            ppfd_10=g.node(vid).Ei10,
            a_n=params.exchange.Na_dict['aN'],
            b_n=params.exchange.Na_dict['bN'],
            a_m=params.exchange.Na_dict['aM'],
            b_m=params.exchange.Na_dict['bM'])

    return gdd_since_budbreak


def remove_stem_geometry(g: MTG):
    geom_prop = g.properties()['geometry']
    to_remove = [i for i in geom_prop if not g.node(i).label.startswith(('L', 'other', 'soil'))]
    [geom_prop.pop(x) for x in to_remove]
    g.properties()['geometry'] = geom_prop
    pass


def set_photosynthetic_capacity(g: MTG, photo_n_params: dict, deactivation_enthalopy: float, leaf_lbl_prefix: str):
    for vid in get_leaves(g=g, leaf_lbl_prefix=leaf_lbl_prefix):
        n = g.node(vid)
        nitrogen_content = n.properties()['Na']
        n.Vcm25 = photo_n_params['Vcm25_N'][0] * nitrogen_content + photo_n_params['Vcm25_N'][1]
        n.Jm25 = photo_n_params['Jm25_N'][0] * nitrogen_content + photo_n_params['Jm25_N'][1]
        n.TPU25 = photo_n_params['TPU25_N'][0] * nitrogen_content + photo_n_params['TPU25_N'][1]
        n.Rd = photo_n_params['Rd_N'][0] * nitrogen_content + photo_n_params['Rd_N'][1]
        n.dHd = deactivation_enthalopy
    pass


def set_collar_water_potential_function(params: Params, **kwargs) -> Callable:
    if params.soil.is_rhyzo_solution:
        def func(transpiration: float, soil_water_potential: float, **kwargs) -> float:
            return soil.calc_collar_water_potential(
                transpiration=transpiration,
                bulk_soil_water_potential=soil_water_potential,
                rhyzosphere_volume=params.soil.rhyzo_volume,
                soil_class=params.soil.soil_class,
                root_radius=params.soil.avg_root_radius,
                root_length=params.soil.root_length)
    else:
        def func(soil_water_potential: float, **kwargs):
            return soil_water_potential
    return func


def init_model(g: MTG, inputs: HydroShootInputs) -> MTG:
    params = inputs.params
    vid_collar = get_mtg_base(g, vtx_label=params.mtg_api.collar_label)
    g.node(g.root).vid_collar = vid_collar
    g.node(g.root).vid_base = vid_collar

    # Add form factors
    if params.simulation.is_energy_budget:
        if inputs.form_factors is not None:
            for ff in ('ff_sky', 'ff_leaves', 'ff_soil'):
                g.properties()[ff] = inputs.form_factors[ff]
        else:
            if not all([s in g.property_names() for s in ('ff_sky', 'ff_leaves', 'ff_soil')]):
                print('Computing form factors...')
                g = set_form_factors_simplified(
                    g=g,
                    pattern=params.irradiance.pattern,
                    infinite=True,
                    leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix,
                    turtle_sectors=params.irradiance.turtle_sectors,
                    icosphere_level=params.irradiance.icosphere_level,
                    unit_scene_length=params.simulation.unit_scene_length)

    # Initialize sap flow to 0
    for vtx_id in traversal.pre_order2(g, g.node(g.root).vid_base):
        g.node(vtx_id).Flux = 0.

    # Add soil surface
    if 'Soil' not in g.properties()['label'].values():
        side_length = inputs.soil_size if inputs.soil_size is not None and inputs.soil_size > 0 else 500.
        g = add_soil_surface_mesh(g=g, side_length=side_length)

    if not inputs.is_ppfd_interception_calculated:
        # Remove undesired geometry for light and energy calculations
        remove_stem_geometry(g)

        # Attach optical properties to MTG elements
        g = set_optical_properties(
            g=g,
            wave_band='SW',
            leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix,
            stem_lbl_prefix=params.mtg_api.stem_lbl_prefix,
            opt_prop=params.irradiance.opt_prop)

    # Calculate leaf Nitrogen per unit surface area according to Prieto et al. (2012)
    if 'Na' not in g.property_names():
        if inputs.leaf_nitrogen is not None:
            g.properties()['Na'] = inputs.leaf_nitrogen
        else:
            print('Computing Nitrogen profile...')
            inputs.gdd_since_budbreak = calc_nitrogen_distribution(
                g=g,
                gdd_since_budbreak=inputs.gdd_since_budbreak,
                weather=inputs.weather,
                params=params)

    set_photosynthetic_capacity(
        g=g,
        photo_n_params=inputs.params.exchange.par_photo_N,
        deactivation_enthalopy=inputs.params.exchange.par_photo['dHd'],
        leaf_lbl_prefix=inputs.params.mtg_api.leaf_lbl_prefix)
    return g


def init_hourly(g: MTG, inputs_hourly: HydroShootHourlyInputs, leaf_ppfd: dict,
                params: Params) -> (MTG, float):
    # Add a date index to g
    g.date = datetime.strftime(inputs_hourly.date, "%Y%m%d%H%M%S")

    # initiate local wind speed
    g.properties()['u'] = set_local_wind_speed(
        g=g, meteo=inputs_hourly.weather, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    # initiate local air temperature
    g.properties()['Tac'] = set_local_air_temperature(
        g=g, meteo=inputs_hourly.weather, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    # Initialize leaf temperature to air temperature
    g = set_leaf_temperature_to_air_temperature(g=g, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    # initiate local leaf-to-air vapor pressure deficit
    g = set_local_vpd(
        g=g, relative_humidity=inputs_hourly.weather['hs'].iloc[0], leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    if leaf_ppfd is not None:
        diffuse_to_total_irradiance_ratio = leaf_ppfd[g.date]['diffuse_to_total_irradiance_ratio']
        g.properties()['Ei'] = leaf_ppfd[g.date]['Ei']
        g.properties()['Eabs'] = leaf_ppfd[g.date]['Eabs']
    else:
        # Compute irradiance distribution over the scene
        caribu_source, diffuse_to_total_irradiance_ratio = irradiance_distribution(
            meteo=inputs_hourly.weather,
            geo_location=params.simulation.geo_location,
            irradiance_unit=params.irradiance.E_type,
            time_zone=params.simulation.tzone,
            turtle_sectors=params.irradiance.turtle_sectors,
            turtle_format=params.irradiance.turtle_format,
            sun2scene=inputs_hourly.sun2scene,
            rotation_angle=params.planting.scene_rotation)

        # Compute irradiance interception and absorbtion
        g, _ = hsCaribu(
            mtg=g,
            unit_scene_length=params.simulation.unit_scene_length,
            source=caribu_source, direct=False,
            infinite=True, nz=50, ds=0.5,
            pattern=params.irradiance.pattern)

    g.properties()['Rg'] = {k: v / (0.48 * 4.6) for k, v in g.properties()['Ei'].items()}

    return g, diffuse_to_total_irradiance_ratio
