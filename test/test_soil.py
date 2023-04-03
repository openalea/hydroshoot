from copy import deepcopy
from itertools import product

from numpy import linspace
from numpy.random import default_rng
from numpy.testing import assert_almost_equal

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

    for soil_class in SOIL_CLASSES:
        sc = soil.SoilTexture(*soil.SOIL_PROPS[soil_class])
        theta_fc = soil.calc_volumetric_water_content_from_water_potential(
            psi=-33 * 10,  # kPa to cm
            theta_res=sc.theta_res,
            theta_sat=sc.theta_sat,
            alpha=sc.alpha,
            n=sc.n)
        assert sc.theta_res < theta_fc < sc.theta_sat


def test_calc_volumetric_water_content_from_water_reservoir():
    assert soil.calc_volumetric_water_content_from_water_reservoir(water_content=1, soil_thickness=1) == 1


def test_calc_soil_water_potential():
    for soil_class in SOIL_CLASSES:
        soil_props = soil.SoilTexture(*soil.SOIL_PROPS[soil_class]).to_dict()
        soil_props.pop('k_sat')

        assert 0 == soil.calc_soil_water_potential(theta=soil_props['theta_sat'], **soil_props)
        assert 1.e-12 > soil.calc_soil_water_potential(theta=soil_props['theta_res'], **soil_props)
        res = []
        for theta in linspace(soil_props['theta_res'], soil_props['theta_sat'], 10):
            res.append(soil.calc_soil_water_potential(theta=theta, **soil_props))
        assert all([x <= y for x, y in zip(res, res[1:])])
    pass


def test_set_root_density_profile():
    for model in ('linear', 'exponential'):
        for a, b in product(default_rng().uniform(-10, 10, 10), default_rng().uniform(0, 10, 10)):
            assert_almost_equal(
                actual=sum(soil.calc_root_density_profile(
                    model=model,
                    a=a,
                    shift=b,
                    soil_layers=list(default_rng().uniform(0, 10, 3)))),
                desired=1,
                decimal=6)

    layers = [1, 1, 1]

    for model in ('linear', 'exponential'):
        res = [soil.calc_root_density_profile(model=model, a=1, shift=v, soil_layers=layers) for v in range(10)]
        assert all([(x[0] - x[-1]) >= (y[0] - y[-1]) for x, y in zip(res, res[1:])])

    for model in ('linear', 'exponential'):
        for a in default_rng().uniform(size=100):
            assert_almost_equal(
                actual=soil.calc_root_density_profile(model=model, a=a, shift=1.e12, soil_layers=layers),
                desired=[1 / 3.] * len(layers),
                decimal=6)

    for a in default_rng().uniform(0, 10, 10):
        assert_almost_equal(
            actual=soil.calc_root_density_profile(model='linear', a=a, shift=0, soil_layers=layers),
            desired=[0.555, 0.333, 0.111],
            decimal=3)


def get_vox():
    class VoxClone(soil.SoilVoxel):
        def __init__(self, **kwargs):
            soil.SoilVoxel.__init__(self, **kwargs)

        def to_dict(self):
            return {k: v for k, v in self.__dict__.items() if not any((k.startswith('__'), callable(v)))}

    soil_class = SOIL_CLASSES[default_rng().integers(low=0, high=len(SOIL_CLASSES) + 1)]
    soil_props = soil.SoilTexture(*soil.SOIL_PROPS[soil_class]).to_dict()
    soil_props.pop('k_sat')
    soil_props['theta_fc'] = soil.calc_volumetric_water_content_from_water_potential(psi=-33 * 10, **soil_props)

    return VoxClone(theta_init=soil_props['theta_fc'], thickness=10, **soil_props)


def test_voxel_water_balance():
    vox = get_vox()
    voxel = deepcopy(vox)
    vox.calc_water_balance(water_input=0, water_uptake=0)
    assert vox.drainage == 0
    assert all([(v1 == v2) for ((k, v1), (_, v2)) in zip(vox.to_dict().items(), voxel.to_dict().items())
                if k != 'drainage'])

    vox.calc_water_balance(water_input=0, water_uptake=10)
    assert vox.drainage == 0
    assert vox.water_content == vox.water_content_res
    assert vox.theta == vox.theta_res

    vox.calc_water_balance(water_input=vox.thickness * (vox.theta_fc - vox.theta_res), water_uptake=0)
    assert vox.drainage == 0
    assert vox.theta == vox.theta_fc
    assert vox.water_content == vox.water_content_fc
    assert_almost_equal(actual=vox.water_potential, desired=-33 * 10, decimal=1)

    theta_prev = vox.theta
    wc_prev = vox.water_content
    vox.calc_water_balance(water_input=1, water_uptake=1)
    assert vox.drainage == 0
    assert vox.theta == theta_prev
    assert vox.water_content == wc_prev

    vox.calc_water_balance(water_input=vox.water_content_fc - vox.water_content + 1, water_uptake=0)
    assert vox.drainage == 1
    assert vox.theta == vox.theta_fc
    assert_almost_equal(actual=vox.water_potential, desired=-33 * 10, decimal=1)

    pass


def get_soil_column(root_density_model: str = 'linear', water_potential_init: float = -1.e12):
    a, shift = (0, 1) if root_density_model == 'cst' else (1, 0)
    model = 'linear' if root_density_model in ('linear', 'cst') else 'exponential'
    return soil.Soil(
        texture_class_name=SOIL_CLASSES[default_rng().integers(low=0, high=len(SOIL_CLASSES) + 1)],
        depth=1550,
        thickness_single_layer=100,
        root_density_model=model,
        root_density_decrease_shape_param=a,
        root_density_at_max_depth=shift,
        water_potential_init=water_potential_init)


def test_soil_water_fills_only_the_first_voxel():
    soil_col = get_soil_column()
    soil_col_init = deepcopy(soil_col)

    soil_col.calc_water_balance(
        water_input=(soil_col.voxels[0].theta_fc - soil_col.voxels[0].theta) * soil_col.voxels[0].thickness,
        transpiration=0)
    assert_almost_equal(actual=soil_col.voxels[0].theta, desired=soil_col.voxels[0].theta_fc, decimal=3)
    assert all([vox.drainage == 0 for vox in soil_col.voxels])
    assert all([(vox1.theta >= vox2.theta) for vox1, vox2 in zip(soil_col.voxels, soil_col.voxels[1:])])
    assert all([vox.theta == vox_init.theta for (vox, vox_init) in zip(soil_col.voxels[1:], soil_col_init.voxels[1:])])
    assert_almost_equal(
        *zip(*[(vox.theta, vox_init.theta) for vox, vox_init in zip(soil_col.voxels[1:], soil_col_init.voxels[1:])]),
        decimal=3)


def test_soil_water_fills_all_voxels_without_deep_drainage():
    soil_col = get_soil_column()
    soil_col.calc_water_balance(
        water_input=sum([((vox.theta_fc - vox.theta) * vox.thickness) for vox in soil_col.voxels]),
        transpiration=0)

    assert_almost_equal(*zip(*[(vox.theta, vox.theta_fc) for vox in soil_col.voxels]), decimal=3)
    assert all([(vox1.drainage >= vox2.drainage) for vox1, vox2 in zip(soil_col.voxels, soil_col.voxels[1:])])
    assert_almost_equal(actual=soil_col.voxels[-1].drainage, desired=0, decimal=3)


def test_soil_deep_drainage_occurs_as_water_input_exceeds_holding_capacity():
    soil_col = get_soil_column()
    soil_col.calc_water_balance(
        water_input=sum([((vox.theta_fc - vox.theta) * vox.thickness) for vox in soil_col.voxels]) + 1,
        transpiration=0)

    assert_almost_equal(actual=soil_col.voxels[-1].drainage, desired=1, decimal=3)


def test_soil_water_is_equally_takenup_by_homogeneous_root_profile():
    soil_col = get_soil_column(water_potential_init=0, root_density_model='cst')
    soil_col.calc_water_balance(
        water_input=0,
        transpiration=sum([((vox.theta - vox.theta_res) * vox.thickness) for vox in soil_col.voxels]))

    assert_almost_equal(*zip(*[(vox.theta, vox.theta_res) for vox in soil_col.voxels]), decimal=3)


def test_soil_water_uptake_follows_root_profile():
    for model in ('linear', 'exponential'):
        soil_col = get_soil_column(water_potential_init=0, root_density_model=model)
        soil_col.calc_water_balance(
            water_input=0,
            transpiration=sum([((vox.theta - vox.theta_res) * vox.thickness) for vox in soil_col.voxels]))

        assert all([(vox1.theta <= vox2.theta) for vox1, vox2 in zip(soil_col.voxels, soil_col.voxels[1:])])
