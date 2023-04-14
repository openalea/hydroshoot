from numpy import linspace, arange

from hydroshoot import soil

SOIL_CLASSES = list(soil.SOIL_PROPS.keys())


def test_calc_soil_water_content_from_water_potential():
    soil_props = dict(
        theta_res=0.045,
        theta_sat=0.430,
        alpha=0.145,
        n=2.68)
    assert soil.calc_volumetric_water_content_from_water_potential(psi=0, **soil_props) == soil_props['theta_sat']
    assert soil.calc_volumetric_water_content_from_water_potential(psi=-1.e12, **soil_props) == soil_props['theta_res']


def test_calc_soil_water_potential():
    for v in soil.SOIL_PROPS.values():
        soil_props = {k: v[i] for i, k in enumerate(('theta_res', 'theta_sat', 'alpha', 'n'))}
        assert 0 == soil.calc_soil_water_potential(theta=soil_props['theta_sat'], **soil_props)
        assert 1.e-12 > soil.calc_soil_water_potential(theta=soil_props['theta_res'], **soil_props)
        res = []
        for theta in linspace(soil_props['theta_res'], soil_props['theta_sat'], 10):
            res.append(soil.calc_soil_water_potential(theta=theta, **soil_props))
        assert all([x <= y for x, y in zip(res, res[1:])])
    pass


def test_soil_conductivity_decreases_as_water_potential_decreases():
    for soil_class in soil.SOIL_PROPS.keys():
        soil_conductivity = [soil.calc_soil_conductivity(psi, soil_class) for psi in arange(0, -3, -0.1)]
        assert all(x >= y for x, y in zip(soil_conductivity, soil_conductivity[1:]))


def test_soil_conductivity_maximum_value_is_greater_for_sand_than_clay():
    soil_conductivity = [soil.calc_soil_conductivity(0., soil_class) for soil_class in ('Clay', 'Sand')]
    assert all(x <= y for x, y in zip(soil_conductivity, soil_conductivity[1:]))


def test_calc_root_soil_resistance_decreases_as_soil_conductivity_increases():
    res = [
        soil.calc_root_soil_resistance(soil_conductivity=v, rhyzosphere_volume=1, root_radius=0.0001, root_length=2000)
        for v in arange(0, 1, 0.01)]
    assert all(x >= y for x, y in zip(res, res[1:]))


def test_calc_root_soil_resistance_decreases_as_rhyzosphere_volume_increases():
    res = [
        soil.calc_root_soil_resistance(soil_conductivity=1, rhyzosphere_volume=v, root_radius=0.0001, root_length=2000)
        for v in arange(.1, 5, 0.01)]
    assert all(x >= y for x, y in zip(res, res[1:]))


def test_calc_root_soil_resistance_decreases_as_root_radius_increases():
    res = [soil.calc_root_soil_resistance(soil_conductivity=1, rhyzosphere_volume=1, root_radius=v, root_length=2000)
           for v in arange(0, 0.001, 0.0001)]
    assert all(x >= y for x, y in zip(res, res[1:]))


def test_calc_root_soil_resistance_decreases_as_root_length_increases():
    res = [soil.calc_root_soil_resistance(soil_conductivity=1, rhyzosphere_volume=1, root_radius=0.0001, root_length=v)
           for v in arange(10, 2000, 10)]
    assert all(x >= y for x, y in zip(res, res[1:]))


def test_calc_collar_water_potential_decreases_as_transpiration_flux_increases():
    root_params = dict(bulk_soil_water_potential=-0.5, rhyzosphere_volume=0.5, root_radius=0.015)
    root_params.update(
        {'root_length': 0.99 * root_params['rhyzosphere_volume'] / (3.14 * (2 * root_params['root_radius']) ** 2)})
    for soil_class in SOIL_CLASSES:
        res = [soil.calc_collar_water_potential(transpiration=v, soil_class=soil_class, **root_params)
               for v in linspace(0, 4 / 3600., 10)]  # max transpiration flux assumed to be 4 L/h
        assert all(x >= y for x, y in zip(res, res[1:]))
