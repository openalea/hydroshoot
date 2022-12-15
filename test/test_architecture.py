from os.path import dirname, join

from hydroshoot import architecture

Source_Dir = join(dirname(__file__), 'data')


def test_cordon():
    """Returns an `openalea.mtg` representing a potted syrah grapevine."""
    g = architecture.vine_mtg(join(Source_Dir, 'digit_vsp.csv'))
    architecture.cordon_vector(g=g)
