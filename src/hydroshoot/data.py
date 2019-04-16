""" A collection of  data objects to be used for testing"""

import json
from pandas import read_csv, DatetimeIndex
from openalea.mtg import traversal
from hydroshoot import architecture
from hydroshoot.data_access import get_path


def potted_syrah():
    """ A mtg representing a potted syrrah grapevine"""
    digit = get_path('grapevine_pot.csv')
    g = architecture.vine_mtg(digit)
    # Local Coordinates Correction
    for v in traversal.iter_mtg2(g, g.root):
        n = g.node(g.Trunk(v, Scale=1)[0])
        theta = 180 if int(n.index()) < 200 else -90 if int(n.index()) < 300 else 0
        architecture.vine_orientation(g, v, theta, local_rotation=True)
    # rotation
    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_orientation(g, v, 90., local_rotation=False)
    for v in traversal.iter_mtg2(g, g.root):
        architecture.vine_phyto_modular(g, v)
        architecture.vine_mtg_properties(g, v)
        architecture.vine_mtg_geometry(g, v)
        architecture.vine_transform(g, v)

    return g


def meteo():
    """A pandas.DataFrame containing meteorological data"""
    meteo_path = get_path('meteo.input')
    df = read_csv(meteo_path, sep=';', decimal='.', header=0)
    df.time = DatetimeIndex(df.time)
    df = df.set_index(df.time)
    return df


def json_parameters():
    params_path = get_path('params.json')
    with open(params_path) as f:
        pars = json.load(f)
    return pars


def reference_time_series_output():
    path = get_path('reference_time_series.output')
    return read_csv(path,sep=';', decimal='.')