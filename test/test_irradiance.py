from hydroshoot.data import potted_syrah, meteo
from hydroshoot.irradiance import irradiance_distribution, hsCaribu, optical_prop, e_conv_PPFD
from numpy.testing import assert_almost_equal


def test_irradiance_distribution():
    # sample values copied from data.json_parameters
    location = (43.61, 3.87, 44.0)
    e_type = 'Rg_Watt/m2'
    conv = e_conv_PPFD(e_type)

    # a cloudy hour
    met = meteo().iloc[[12], :]
    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs == 1
    # would expect 46 as no direct light is estimated
    assert len(sources) == 47
    nrj, pos = zip(*sources)
    assert_almost_equal(met.Rg.sum() * conv, sum(nrj), 2)

    # Same our but with 16 directions
    sources, rdrs = irradiance_distribution(met, location, e_type, turtle_sectors='16')
    assert len(sources) == 17
    nrj, pos = zip(*sources)
    assert_almost_equal(met.Rg.sum() * conv, sum(nrj), 2)

    # a sunny hour
    met = meteo().iloc[[60], :]
    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs < 0.3
    assert len(sources) == 47
    nrj, pos = zip(*sources)
    assert_almost_equal(met.Rg.sum() * conv, sum(nrj), 2)

    #a night hour
    met = meteo().iloc[[1], :]
    sources, rdrs = irradiance_distribution(met, location, e_type)
    assert rdrs == 1
    # would expect No source or 46
    assert len(sources) == 47
    nrj, pos = zip(*sources)
    assert_almost_equal(sum(nrj), 0)

    # a particularly overcast day
    day_met = meteo().iloc[:24,:]
    sources, rdrs = irradiance_distribution(day_met, location, e_type)
    assert rdrs > 0.9
    # not optimals  as almost no direct light is predicted for this day. Expect 48 instead
    assert len(sources) == 70
    nrj, pos = zip(*sources)
    assert_almost_equal(day_met.Rg.sum() * conv, sum(nrj),2)

    # Same day but with 16 directions
    sources, rdrs = irradiance_distribution(day_met, location, e_type, turtle_sectors='16')
    assert len(sources) == 40
    nrj, pos = zip(*sources)
    assert_almost_equal(day_met.Rg.sum() * conv, sum(nrj), 2)

    # a much clearer day
    day_met = meteo().iloc[48:72,:]
    sources, rdrs = irradiance_distribution(day_met, location, e_type)
    assert rdrs < 0.5
    # not optimals
    assert len(sources) == 70
    nrj, pos = zip(*sources)
    assert_almost_equal(day_met.Rg.sum() * conv,sum(nrj),0)


def test_hsCaribu():
    g = potted_syrah()
    # sample values copied from data.json_parameters
    location = (43.61, 3.87, 44.0)
    e_type = 'Rg_Watt/m2'
    unit_scene_length = 'cm'

    # reporduce preprocessing as done  in model.py
    # Suppression of undesired geometry for light and energy calculations
    geom_prop = g.properties()['geometry']
    vidkeys = []
    for vid in g.properties()['geometry']:
        n = g.node(vid)
        if not n.label.startswith(('L', 'other', 'soil')):
            vidkeys.append(vid)
    [geom_prop.pop(x) for x in vidkeys]
    g.properties()['geometry'] = geom_prop
    # Attaching optical properties to MTG elements
    g = optical_prop(g)


    # force sources
    day_met = None
    local_date = None
    sources = [(1,(0,0,-1))]
    # pattern = False (default) raises a bug
    gg, cs = hsCaribu(g, day_met, local_date, location, e_type, unit_scene_length, source=sources)
