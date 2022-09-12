from datetime import datetime, timedelta

import openalea.mtg.traversal as traversal
from openalea.mtg.mtg import MTG
from pandas import DataFrame
from pandas import date_range

from hydroshoot.architecture import get_mtg_base, add_rhyzosphere_concentric_cylinders, add_soil_surface_mesh
from hydroshoot.energy import set_form_factors_simplified, set_wind_speed
from hydroshoot.exchange import leaf_Na
from hydroshoot.io import HydroShootInputs, HydroShootHourlyInputs
from hydroshoot.irradiance import irradiance_distribution, hsCaribu
from hydroshoot.irradiance import set_optical_properties
from hydroshoot.params import Params
from hydroshoot.preprocess import calc_gdd_since_budbreak


def calc_nitrogen_distribution(g: MTG, gdd_since_budbreak: float, weather: DataFrame, params: Params) -> (MTG, float):
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
        rotation_angle=params.irradiance.scene_rotation)

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
    for vid in g.VtxList(Scale=3):
        if g.node(vid).label.startswith(params.mtg_api.leaf_lbl_prefix):
            g.node(vid).Na = leaf_Na(
                age_gdd=gdd_since_budbreak,
                ppfd_10=g.node(vid).Ei10,
                a_n=params.exchange.Na_dict['aN'],
                b_n=params.exchange.Na_dict['bN'],
                a_m=params.exchange.Na_dict['aM'],
                b_m=params.exchange.Na_dict['bM'])

    return g, gdd_since_budbreak


def remove_stem_geometry(g: MTG) -> MTG:
    geom_prop = g.properties()['geometry']
    to_remove = [i for i in geom_prop if not g.node(i).label.startswith(('L', 'other', 'soil'))]
    [geom_prop.pop(x) for x in to_remove]
    g.properties()['geometry'] = geom_prop

    return g


def init_model(g: MTG, inputs: HydroShootInputs) -> MTG:
    params = inputs.params
    g.node(g.root).vid_collar = get_mtg_base(g, vtx_label=params.mtg_api.collar_label)
    g.node(g.root).vid_base = g.node(g.root).vid_collar

    # Add rhyzosphere concentric cylinders
    if params.soil.rhyzo_solution:
        if not any(item.startswith('rhyzo') for item in g.property('label').values()):
            g = add_rhyzosphere_concentric_cylinders(
                g=g,
                cylinders_radii=params.soil.rhyzo_radii,
                soil_dimensions=params.soil.soil_dimensions,
                soil_class=params.soil.soil_class,
                vtx_label=params.mtg_api.collar_label,
                length_conv=1. / params.simulation.conv_to_meter)

    # Add form factors
    if params.simulation.energy_budget:
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

    # Remove undesired geometry for light and energy calculations
    g = remove_stem_geometry(g)

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
            g, inputs.gdd_since_budbreak = calc_nitrogen_distribution(
                g=g,
                gdd_since_budbreak=inputs.gdd_since_budbreak,
                weather=inputs.weather,
                params=params)

    return g


def init_hourly(g: MTG, inputs_hourly: HydroShootHourlyInputs, leaf_absorbed_ppfd: dict,
                params: Params) -> (MTG, float):
    # Add a date index to g
    g.date = datetime.strftime(inputs_hourly.date, "%Y%m%d%H%M%S")

    # initiate wind speed
    g.properties()['u'] = set_wind_speed(
        g=g, meteo=inputs_hourly.weather, leaf_lbl_prefix=params.mtg_api.leaf_lbl_prefix)

    if leaf_absorbed_ppfd is not None:
        g.properties()['Ei'], diffuse_to_total_irradiance_ratio = leaf_absorbed_ppfd[inputs_hourly.date]
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
            rotation_angle=params.irradiance.scene_rotation)

        # Compute irradiance interception and absorbtion
        g, _ = hsCaribu(
            mtg=g,
            unit_scene_length=params.simulation.unit_scene_length,
            source=caribu_source, direct=False,
            infinite=True, nz=50, ds=0.5,
            pattern=params.irradiance.pattern)

    return g, diffuse_to_total_irradiance_ratio
