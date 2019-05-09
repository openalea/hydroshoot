from non_regression_data import potted_syrah, meteo
from hydroshoot.energy import form_factors_simplified, leaf_temperature, forced_soil_temperature
from numpy.testing import assert_almost_equal
import openalea.plantgl.all as pgl
import hydroshoot.energy as energy
import warnings


def test_pgl_scene():
    g = potted_syrah()
    bc = pgl.BBoxComputer(pgl.Tesselator())

    s = energy.pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(),bbox.getZMax()
    assert zmax > zmin > 0

    s = energy.pgl_scene(g, flip=True)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(),bbox.getZMax()
    assert 0 > zmax > zmin

    # check that original scene is still z >0
    s = energy.pgl_scene(g)
    bc.process(s)
    bbox = bc.boundingbox
    zmin, zmax = bbox.getZMin(), bbox.getZMax()
    assert zmax > zmin > 0


def test_get_leaves():
    g = potted_syrah()
    leaves = energy.get_leaves(g, leaf_lbl_prefix='L')
    assert len(leaves) == 46


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


def test_heat_boundary_layer_conductance():
    g = potted_syrah()
    met = meteo().iloc[[12], :]
    gbH = energy.heat_boundary_layer_conductance(g, met)
    assert len(gbH) == 46
    assert_almost_equal(sum(gbH.values()) / len(gbH), 47, 0)


def test_forced_soil_temperature():
    met = meteo().iloc[[12], :]
    tsoil = forced_soil_temperature(met)
    assert tsoil == met.Tac[0] + 20


def test_leaf_temperature():
    g = potted_syrah()
    met = meteo().iloc[[12], :]
    tsoil = 20
    tsky = 2

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('default')
        tleaf, it = leaf_temperature(g, met, tsoil, tsky)
        assert len(w) == 1
        assert len(tleaf) == 46
        for vid in tleaf:
            assert tleaf[vid] == met.Tac[0]

    g.properties()['gbH'] = energy.heat_boundary_layer_conductance(g, met)
    tleaf, it = leaf_temperature(g, met, tsoil, tsky, k_sky=0.5, k_soil=0.5, k_leaves=1, Ei=0, E=0)
    assert len(tleaf) == 46
    for vid in tleaf:
        assert tleaf[vid] != met.Tac[0]