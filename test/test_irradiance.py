from numpy.testing import assert_almost_equal

from openalea.hydroshoot.irradiance import irradiance_distribution, hsCaribu, set_optical_properties, e_conv_PPFD
from test.non_regression_data import potted_syrah, meteo


def test_irradiance_distribution():
    # sample values copied from data.json_parameters
    location = (43.61, 3.87, 44.0)
    e_type = 'Rg_Watt/m2'
    conv = e_conv_PPFD(e_type)

    # a cloudy hour
    met = meteo().iloc[[12], :]
    rg = met.Rg.sum() * conv

    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs == 1
    assert len(sources) == 46
    nrj, pos = zip(*sources)
    assert abs(rg - sum(nrj)) / rg < 0.001

    # Same our but with 16 directions
    sources, rdrs = irradiance_distribution(met, location, e_type, turtle_sectors='16')
    assert len(sources) == 16
    nrj, pos = zip(*sources)
    assert abs(rg - sum(nrj)) / rg < 0.001

    # a sunny hour
    met = meteo().iloc[[60], :]
    rg = met.Rg.sum() * conv

    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs < 0.3
    assert len(sources) == 47
    nrj, pos = zip(*sources)
    assert abs(rg - sum(nrj)) / rg < 0.001

    # a night hour
    met = meteo().iloc[[1], :]
    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs == 1
    assert len(sources) == 1
    nrj, pos = zip(*sources)
    assert_almost_equal(sum(nrj), 0)

    # a particularly overcast day
    day_met = meteo().iloc[:24, :]
    rg_total = day_met.Rg.sum() * conv

    sources, rdrs = irradiance_distribution(day_met, location, e_type)
    assert rdrs > 0.9
    assert len(sources) == 47
    nrj, pos = zip(*sources)
    assert abs(rg_total - sum(nrj)) / rg_total < 0.001

    # Same day but with 16 directions
    sources, rdrs = irradiance_distribution(day_met, location, e_type, turtle_sectors='16')
    assert len(sources) == 17
    nrj, pos = zip(*sources)
    assert abs(rg_total - sum(nrj)) / rg_total < 0.001

    # a much clearer day
    day_met = meteo().iloc[48:72, :]
    rg_total = day_met.Rg.sum() * conv
    sources, rdrs = irradiance_distribution(day_met, location, e_type)
    assert rdrs < 0.5
    assert len(sources) == 58
    nrj, pos = zip(*sources)
    assert abs(rg_total - sum(nrj)) / rg_total < 0.001


def test_hsCaribu():
    args = dict(
        g=potted_syrah(),
        wave_band='SW',
        leaf_lbl_prefix='L',
        stem_lbl_prefix=('in', 'Pet', 'cx'),
        opt_prop={'SW': {'leaf': (0.06, 0.07), 'stem': (0.13,), 'other': (0.65, 0.0)},
                  'LW': {'leaf': (0.04, 0.07), 'stem': (0.13,), 'other': (0.65, 0.0)}})
    # sample values copied from data.json_parameters
    unit_scene_length = 'cm'
    # Attaching optical properties to MTG elements
    g = set_optical_properties(**args)

    # simple run
    assert 'Ei' not in g.property_names()
    assert 'Eabs' not in g.property_names()
    g, cs = hsCaribu(g, unit_scene_length)
    assert 'Ei' in g.property_names()
    assert 'Eabs' in g.property_names()
    assert len(g.property('Ei')) == len(cs.scene)
    assert len(g.property('Eabs')) == len(cs.scene)
    assert sum(g.property('Ei').values()) > 0

    # simple run during night
    g, cs = hsCaribu(g, unit_scene_length, source=[(0, (0, -1, 0))])
    assert sum(g.property('Ei').values()) == 0

    # reproduce scene filtering as done in model.py
    args.update({'g': potted_syrah()})
    g = set_optical_properties(**args)
    ng = len(g.property('geometry'))
    label = g.property('label')
    g.properties()['radiative_geometry'] = {k: v for k, v in g.property('geometry').items() if
                                            label[k].startswith(('L', 'other', 'soil'))}
    assert len(g.property('radiative_geometry')) < ng
    g, cs = hsCaribu(g, unit_scene_length, geometry='radiative_geometry')
    assert len(g.property('Ei')) == len(cs.scene)
    assert len(g.property('Ei')) == len(g.property('radiative_geometry'))
    assert len(g.property('geometry')) == ng
    # non regression test
    ei_sum = sum(g.property('Ei').values())
    assert_almost_equal(ei_sum, 14.83, 2)

    # test consider option
    args.update({'g': potted_syrah()})
    g = set_optical_properties(**args)
    ng = len(g.property('geometry'))
    label = g.property('label')
    consider = [vid for vid in g.property('geometry') if label[vid].startswith(('L', 'other', 'soil'))]
    assert len(consider) < ng
    g, cs = hsCaribu(g, unit_scene_length, consider=consider)
    assert len(g.property('Ei')) == len(cs.scene)
    assert len(g.property('Ei')) == len(consider)
    assert len(g.property('geometry')) == ng
    # non regression test
    ei_sum = sum(g.property('Ei').values())
    assert_almost_equal(ei_sum, 14.83, 2)
