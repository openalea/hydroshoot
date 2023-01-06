# -*- coding: utf-8 -*-
"""Irradiance module of HydroShoot: an interface of the Caribu module.
This module computes leaf (and eventually other elements) interception and absorption of incident irradiance given
plant shoot geometry.

TODO: plug to the standard interface of Caribu module.
"""

import alinea.astk.icosphere as ico
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import turtle, Gensun, GetLightsSun
from alinea.caribu.sky_tools.spitters_horaire import RdRsH
from numpy import array, deg2rad
from openalea.plantgl.all import Translated, Sphere, Shape, Material, Color3, Viewer
from pandas import date_range
from pvlib.solarposition import ephemeris
from pytz import timezone

from hydroshoot.architecture import vector_rotation


def local2solar(local_time, latitude, longitude, time_zone, temperature=25.):
    """Calculates UTC time and Solar time in decimal hours (solar noon is 12.00), based on
    :func:`pvlib.solarposition.ephemeris`.

    Args:
        local_time (datetime): local time
        latitude (float): [°] latitude at which local time is given
        longitude (float): [°] longitude at which local time is given
        time_zone (str): a 'pytz.timezone' (e.g. 'Europe/Paris')
        temperature (float): [°C] air temperature

    Returns:
        (float): [decimal hours] solar time

    """

    time_zone = timezone(time_zone)
    local_time = time_zone.localize(local_time)
    #    date_utc = local_time.astimezone(utc)
    datet = date_range(local_time, local_time)
    solar_time = ephemeris(datet, latitude, longitude, temperature).solar_time.values[0]

    return solar_time


def set_optical_properties(g, leaf_lbl_prefix, stem_lbl_prefix, wave_band, opt_prop):
    """Attaches optical properties to elements of an MTG. An element may be a `leaf`, a `stem`, or `other` when scene
    background elements (e.g. oculting objects) are attached to the MTG.

    Args:
        g (MTG): plant Multiscale Tree Graph
        leaf_lbl_prefix (str): the prefix of the leaf label (e.g. 'L')
        stem_lbl_prefix (str or tuple of str): the prefix of the stem label (e.g. ('in', 'Pet', 'cx'))
        wave_band (str): wave band name, default values are 'SW' for *short wave band* or 'LW' for *long wave band*
        opt_prop: dictionary of the optical properties of mtg elements,
            given as a {band_name: material} dictionary of tuples
            e.g. : {'SW': {'leaf': (0.06, 0.07), 'stem': (0.13,), 'other': (0.65, 0.0)},
                    'LW': {'leaf': (0.04, 0.07), 'stem': (0.13,), 'other': (0.65, 0.0)}}
            (See :func:`CaribuScene.__init__` for more information)

    Returns:
        (MTG): mtg object

    """
    leaf_material, stem_material = [opt_prop[wave_band][ikey] for ikey in ('leaf', 'stem')]
    other_material = opt_prop[wave_band]['other']

    geom = g.property('geometry')

    for vid in geom:
        n = g.node(vid)
        if n.label.startswith(leaf_lbl_prefix):
            n.opticals = leaf_material
        elif n.label.startswith(stem_lbl_prefix):
            n.opticals = stem_material
        elif n.label.startswith('other'):
            n.opticals = other_material

    return g


def e_conv_Wm2(irradiance_unit):
    """Conversion factor for irradiance from the unit specified by `irradiance_unit` to [W m-2]

    Args:
        irradiance_unit (str): unit of the irradiance flux density,
            one of ('Rg_Watt/m2', 'RgPAR_Watt/m2', 'PPFD_umol/m2/s')

    Returns:
        (float): conversion factor

    """

    if irradiance_unit == 'Rg_Watt/m2':
        return 1.
    elif irradiance_unit == 'RgPAR_Watt/m2':
        return 1. / 0.48
    elif irradiance_unit == 'PPFD_umol/m2/s':
        return 1. / (0.48 * 4.6)
    else:
        raise TypeError("E_type must be one of the following 'Rg_Watt/m2', 'RgPAR_Watt/m2' or'PPFD_umol/m2/s'.")


def e_conv_PPFD(irradiance_unit):
    """Conversion factor for irradiance from the unit specified by `irradiance_unit` to [umol/m2/s]

    Args:
        irradiance_unit (str): unit of the irradiance flux density,
            one of ('Rg_Watt/m2', 'RgPAR_Watt/m2', 'PPFD_umol/m2/s')

    Returns:
        (float): conversion factor

    """

    if irradiance_unit == 'Rg_Watt/m2':
        return 0.48 * 4.6
    elif irradiance_unit == 'RgPAR_Watt/m2':
        return 4.6
    elif irradiance_unit == 'PPFD_umol/m2/s':
        return 1.
    else:
        raise TypeError("E_type must be one of the following 'Rg_Watt/m2', 'RgPAR_Watt/m2' or'PPFD_umol/m2/s'.")


def irradiance_distribution(meteo, geo_location, irradiance_unit,
                            time_zone='Europe/Paris', turtle_sectors='46', turtle_format='uoc',
                            sun2scene=None, rotation_angle=0., icosphere_level=None):
    """Calculates irradiance distribution over a semi-hemisphere surrounding the plant [umol m-2 s-1].

    Args:
        meteo (DataFrame): meteo data having the following columns:
            time (datetime): UTC time
            Tac (float): [°C] air temperature
            hs (float): (%) air relative humidity
            Rg or PPFD (float): respectively global [W m-2] or photosynthetic photon flux density [umol m-2 s-1]
            u (float): [m s-1] wind speed
            Ca (float): [ppm] CO2 concentration in the air
            Pa (float): [kPa] atmospheric pressure
        geo_location: tuple of (latitude [°], longitude [°], elevation [°])
        irradiance_unit (str): unit of the irradiance flux density,
            one of ('Rg_Watt/m2', 'RgPAR_Watt/m2', 'PPFD_umol/m2/s')
        time_zone (str): a 'pytz.timezone' (e.g. 'Europe/Paris')
        turtle_sectors (str): number of turtle sectors (see :func:`turtle` from `sky_tools` package)
        turtle_format (str): format irradiance distribution, could be 'soc', or 'uoc'
            (see :func:`turtle` from `sky_tools` package for details)
        sun2scene (pgl.scene): if provided, a sun object (sphere) is added to it
        rotation_angle (float): [°] counter clockwise angle between the default X-axis direction (South) and real
            direction of X-axis
        icosphere_level (int): the level of refinement of the dual icosphere
            (see :func:`alinea.astk.icosphere.turtle_dome` for details)

    Returns:
        [umol m-2 s-1] tuple of tuples, cumulative irradiance flux densities distributed across the semi-hemisphere
            surrounding the plant
        (float) [-] diffuse-to-total irradiance ratio

    Notes:
        meteo data can consist of only one line (single event) or multiple lines.
            In the latter case, this function returns accumulated irradiance throughtout the entire periode with sun
            positions corresponding to each time step.


    TODO: replace by the icosphere procedure

    """
    diffuse_ratio = []
    nrj_sum = 0
    for date, weather in meteo.iterrows():

        energy = meteo.loc[date, 'PPFD' if irradiance_unit.split('_')[0] == 'PPFD' else 'Rg']

        # First check: convert irradiance to W m-2 (Spitters method always gets energy flux as Rg Watt m-2)
        corr = e_conv_Wm2(irradiance_unit)
        energy = energy * corr

        # Second check: Convert to solar datetime
        doy = date.dayofyear
        latitude, longitude, elevation = [geo_location[x] for x in range(3)]
        temperature = meteo.Tac.values[0]
        solar_time = local2solar(date, latitude, longitude, time_zone, temperature)
        # doy_utc = date_utc.timetuple().tm_yday
        # hour_utc = date_utc.hour + date_utc.minute / 60.

        diffuse_ratio_hourly = RdRsH(Rg=energy, DOY=doy, heureTU=solar_time, latitude=latitude)
        diffuse_ratio.append(diffuse_ratio_hourly * energy)
        nrj_sum += energy

        # Third and final check: it is always desirable to get energy as PPFD
        energy = energy * (0.48 * 4.6)

        irradiance_diff = diffuse_ratio_hourly * energy
        irradiance_dir = (1 - diffuse_ratio_hourly) * energy

        # diffuse irradiance
        if not icosphere_level:
            energy2, emission, direction, elevation, azimuth = turtle.turtle(sectors=turtle_sectors,
                                                                             format=turtle_format,
                                                                             energy=irradiance_diff)
        else:
            vert, fac = ico.turtle_dome(icosphere_level)
            direction = ico.sample_faces(vert, fac, iter=None, spheric=False).values()
            direction = [idirect[0] for idirect in direction]
            direction = map(lambda x: tuple(list(x[:2]) + [-x[2]]), direction)

        sky = list(zip(energy2, direction))

        # direct irradiance
        sun = Gensun.Gensun()(Rsun=irradiance_dir, DOY=doy, heureTU=solar_time, lat=latitude)
        sun = GetLightsSun.GetLightsSun(sun)
        sun_data = [(float(sun.split()[0]), (float(sun.split()[1]), float(sun.split()[2]), float(sun.split()[3])))]

        # diffuse irradiance (distributed over a dome) + direct irradiance (localized as point source(s))
        source = sky.__add__(sun_data)
        source = [list(isource) for isource in source]

        try:
            for j in range(len(source) - 1):
                source_cum[j][0] += source[j][0]
            source_cum.append(source[-1])
        except NameError:
            source_cum = []
            for isource in source:
                source_cum.append([isource[0], isource[1]])

        if date == meteo.index[-1]:
            source_cum = [tuple(isource) for isource in source_cum]

    # Rotate irradiance sources to cope with plant row orientation
    if rotation_angle != 0.:
        v_energy = [vec[0] for vec in source_cum]
        v_coord = [tuple(vector_rotation(vec[1], (0., 0., 1.), deg2rad(rotation_angle))) for vec in source_cum]
        source_cum = zip(v_energy, v_coord)

    # Add Sun to an existing pgl.scene
    if sun2scene is not None:
        xSun, ySun, zSun = -500. * array([source_cum[-1][1][i] for i in range(3)])
        if zSun >= 0:
            ss = Translated(xSun, ySun, zSun, Sphere(20))
            sun = Shape(ss, Material('yellow', Color3(255, 255, 0)))
            sun2scene.add(sun)
        Viewer.display(sun2scene)

    # Diffuse_ratio mean
    if nrj_sum > 0:
        diffuse_ratio = sum(diffuse_ratio) / nrj_sum
    else:
        diffuse_ratio = 1

    # Filter black sources up to the penultimate (hsCaribu expect at least one source)
    source_cum = [v for v in source_cum if v[0] > 0]
    if len(source_cum) == 0:  # night
        source_cum = [(0, (0, 0, -1))]

    return source_cum, diffuse_ratio


def hsCaribu(mtg, unit_scene_length, geometry='geometry', opticals='opticals', consider=None,
             source=None, direct=True,
             infinite=False,
             nz=50, ds=0.5, pattern=None, soil_reflectance=0.15):
    """Calculates intercepted and absorbed irradiance flux densities by the plant canopy.

    Args:
        mtg (MTG): plant Multiscale Tree Graph
        unit_scene_length (str): the unit of length used for scene coordinate and for pattern
            (should be one of `CaribuScene.units` default)
        geometry (str): the name of the property to use for computing scene geometry from the mtg
        opticals (str): the name of the property to use for caribu optical properties
        consider (list(int)): vertices to be considered for the computation, if None (default) all vertices with a
            geometry are considered
        source (list): a tuple of tuples, giving energy unit and sky coordinates, if None, a unit zenital source is
            used
        direct: see :func:`runCaribu` from `CaribuScene` package
        infinite: see :func:`runCaribu` from `CaribuScene` package
        nz: see :func:`runCaribu` from `CaribuScene` package
        ds: see :func:`runCaribu` from `CaribuScene` package
        pattern: see :func:`runCaribu` from `CaribuScene` package
        soil_reflectance (float): [-] the reflectance of the soil (between 0 and 1)

    Returns:
        mtg object with the incident irradiance (`Ei`) and absorbed irradiance (`Eabs`), both in [umol m-2 s-1],
            attached to mtg vertices as properties.

    Notes:
        `Ei` and `Eabs` units are returned in [umol m-2 s-1] **REGARDLESS** of the `unit_scence_length` type.

    """
    assert geometry in mtg.property_names()
    assert opticals in mtg.property_names()

    # Currently used as a dummy variable
    wave_band = 'SW'

    if source is None:
        source = [(1, (0, 0, -1))]

    # Hack : would be much better to have geometry as a caribuscene args directly
    geom0 = mtg.property(geometry)
    geometry0 = {}
    remove_geometry = 'geometry' not in mtg.property_names()

    if consider is not None:
        geom0 = {k: v for k, v in mtg.property(geometry).items()}
        mtg.properties()[geometry] = {k: v for k, v in mtg.property(geometry).items() if k in consider}
    if geometry != 'geometry' and 'geometry' in mtg.property_names():
        geometry0 = {k: v for k, v in mtg.property('geometry').items()}

    mtg.properties()['geometry'] = mtg.property(geometry)
    caribu_scene = None

    # Use a try/except/finally block to allow restoring geometry if block fails
    try:
        # Run caribu only if irradiance is greater that 0
        if sum([x[0] for x in source]) == 0.:
            mtg.properties()['Ei'] = {k: 0. for k in mtg.property('geometry')}
            mtg.properties()['Eabs'] = {k: 0. for k in mtg.property('geometry')}
        else:
            # Attaching optical properties to each organ of the plant mock-up
            opts = {wave_band: mtg.property(opticals)}

            # Setup CaribuScene
            caribu_scene = CaribuScene(mtg, light=source, opt=opts,
                                       soil_reflectance={wave_band: soil_reflectance},
                                       scene_unit=unit_scene_length,
                                       pattern=pattern)

            # Run caribu
            raw, aggregated = caribu_scene.run(direct=direct, infinite=infinite, d_sphere=ds, layers=nz,
                                               split_face=False)

            # Attaching output to MTG
            mtg.properties()['Ei'] = aggregated[wave_band]['Ei']
            mtg.properties()['Eabs'] = aggregated[wave_band]['Eabs']
    except:
        pass
    finally:
        mtg.properties()[geometry] = geom0
        if remove_geometry:
            mtg.remove_property('geometry')
        else:
            if len(geometry0) > 0:
                mtg.properties()['geometry'] = geometry0

    return mtg, caribu_scene
