""" A collection of  data objects to be used for testing"""

from os.path import dirname, join
from json import load
from pandas import read_csv, DatetimeIndex

from openalea.mtg import traversal
from hydroshoot import architecture

sources_dir = join(dirname(__file__), 'data')


def potted_syrah():
    """Returns an `openalea.mtg` representing a potted syrah grapevine."""
    digit = join(sources_dir, 'grapevine_pot.csv')
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
    """Returns a `pandas.DataFrame` containing meteorological data."""
    meteo_path = join(sources_dir, 'meteo.input')
    df = read_csv(meteo_path, sep=';', decimal='.', header=0)
    df.time = DatetimeIndex(df.time)
    df = df.set_index(df.time)
    return df


def json_parameters():
    """Returns `Dict` of hydroshoot parameters."""
    params_path = join(sources_dir, 'params.json')
    with open(params_path) as f:
        pars = load(f)
    return pars


def reference_time_series_output():
    """Returns a `pandas.DataFrame` containing reference hydroshoot time-series output."""
    path = join(sources_dir, 'reference_time_series.output')
    return read_csv(path, sep=';', decimal='.')
