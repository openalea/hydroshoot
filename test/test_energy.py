from non_regression_data import potted_syrah, meteo
from hydroshoot.energy import pgl_scene, form_factors_simplified, leaf_temperature_init, leaf_temperature, forced_soil_temperatue
from numpy.testing import assert_almost_equal
import openalea.plantgl.all as pgl


def test_pgl_scene():
    g = potted_syrah()
    bc = pgl.BBoxComputer(pgl.Tesselator())

    s = pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(),bbox.getZMax()
    assert zmax > zmin > 0

    s = pgl_scene(g, flip=True)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(),bbox.getZMax()
    assert 0 > zmax > zmin

    # check that original scene is still z >0
    s = pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(), bbox.getZMax()
    assert zmax > zmin > 0


def test_form_factors_simplified():
    g = potted_syrah()
    assert 'k_soil' not in g.property_names()
    assert 'k_sky' not in g.property_names()
    assert 'k_leaves' not in g.property_names()
    g = form_factors_simplified(g, icosphere_level=0)
    assert 'k_soil' in g.property_names()
    assert 'k_sky' in g.property_names()
    assert 'k_leaves' in g.property_names()
    # non regression test
    assert_almost_equal(sum(g.property('k_leaves').values()), 147.7, 1)


def test_forced_soil_temperature():
    met = meteo().iloc[[12], :]
    tsoil = forced_soil_temperatue(met)
    assert tsoil == met.Tac[0] + 20

def test_leaf_temperature():
    g = potted_syrah()
    met = meteo().iloc[[12], :]
    tsoil = 20
    tsky = 2
    g = leaf_temperature_init(g)
    it = leaf_temperature(g, met, tsoil, tsky)