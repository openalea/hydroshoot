from os.path import dirname, join

from openalea.hydroshoot import architecture

Source_Dir = join(dirname(__file__), 'data')


def test_cordon():
    g = architecture.vine_mtg(join(Source_Dir, 'digit_vsp.csv'))
    architecture.cordon_vector(g=g)


def test_multiple_plants_digit():
    try:
        architecture.vine_mtg(file_path=join(Source_Dir, 'digit_twin_pots.csv'))
    except ValueError:
        pass
