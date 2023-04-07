import openalea.plantgl.all as pgl
from numpy.testing import assert_almost_equal

import hydroshoot.energy as energy
from hydroshoot.architecture import get_leaves
from hydroshoot.energy import set_form_factors_simplified, calc_leaf_temperature, force_soil_temperature
from test.non_regression_data import potted_syrah, meteo


def test_pgl_scene():
    g = potted_syrah()
    bc = pgl.BBoxComputer(pgl.Tesselator())

    s = energy.pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(), bbox.getZMax()
    assert zmax > zmin > 0

    s = energy.pgl_scene(g, flip=True)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(), bbox.getZMax()
    assert 0 > zmax > zmin

    # check that original scene is still z >0
    s = energy.pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(), bbox.getZMax()
    assert zmax > zmin > 0


def test_get_leaves():
    g = potted_syrah()
    leaves = get_leaves(g, leaf_lbl_prefix='L')
    assert len(leaves) == 46


def test_form_factors_simplified():
    g = potted_syrah()
    set_form_factors_simplified(g, icosphere_level=0)
    # non regression test
    assert_almost_equal(sum(g.property('ff_leaves').values()), 148.6, 1)


def test_forced_soil_temperature():
    met = meteo().iloc[[12], :]
    tsoil = force_soil_temperature(met)
    assert tsoil == met.Tac[0] + 20


def test_leaf_temperature():
    g = potted_syrah()
    met = meteo().iloc[[12], :]
    tsoil = 20
    tsky = 2

    for vid in g.properties()['geometry'].keys():
        node = g.node(vid)
        node.Ei = 0
        node.Rg = 0
        node.ff_sky = 0.3
        node.ff_leaves = 0.3
        node.ff_soil = 0.4
        node.gbH = 1
        node.E = 0.
        node.Tlc = met.Tac.values[0]

    tleaf, it = calc_leaf_temperature(g, met, tsoil, tsky, get_leaves(g=g, leaf_lbl_prefix='L'))
    assert len(tleaf) == 46
    first = list(tleaf.keys())[0]
    for vid in tleaf:
        assert tleaf[vid] == tleaf[first]
        assert tleaf[vid] != met.Tac[0]
