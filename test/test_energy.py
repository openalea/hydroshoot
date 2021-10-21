from test.non_regression_data import potted_syrah, meteo
from hydroshoot.energy import form_factors_simplified, leaf_temperature, forced_soil_temperature
from numpy.testing import assert_almost_equal
import openalea.plantgl.all as pgl
import hydroshoot.energy as energy


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
    k_soil, k_sky, k_leaves = form_factors_simplified(g, icosphere_level=0)
    # non regression test
    assert_almost_equal(sum(k_leaves.values()), 147.7, 1)


def test_heat_boundary_layer_conductance():
    g = potted_syrah()
    met = meteo().iloc[[12], :]
    l = energy.get_leaves_length(g)
    gbH = energy.heat_boundary_layer_conductance(l, met.u[0])
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

    tleaf, it = leaf_temperature(g, met, tsoil, tsky)
    assert len(tleaf) == 46
    first = list(tleaf.keys())[0]
    for vid in tleaf:
        assert tleaf[vid] == tleaf[first]
        assert tleaf[vid] != met.Tac[0]

    l = energy.get_leaves_length(g)
    u = energy.leaf_wind_as_air_wind(g, met)
    gbH = energy.heat_boundary_layer_conductance(l, u)
    tleaf, it = leaf_temperature(g, met, tsoil, tsky, gbh=gbH)
    first = list(tleaf.keys())[0]
    assert len(tleaf) == 46
    for vid in tleaf:
        assert tleaf[vid] != met.Tac[0]
        if vid != first:
            assert tleaf[vid] != tleaf[first]