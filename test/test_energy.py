from non_regression_data import potted_syrah, meteo
from hydroshoot.energy import pgl_scene, form_factors_simplified
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
    g = form_factors_simplified(g)
    assert 'k_soil' in g.property_names()
    assert 'k_sky' in g.property_names()
    assert 'k_leaves' in g.property_names()
    # non regression test
    assert_almost_equal(sum(g.property('k_leaves').values()), 148.4, 1)

