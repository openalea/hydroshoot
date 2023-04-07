from datetime import datetime
from pathlib import Path

from openalea.mtg.mtg import MTG
from openalea.plantgl.all import Scene
from pandas import DataFrame, read_csv, DatetimeIndex

from hydroshoot.architecture import get_leaves
from hydroshoot.display import visu
from hydroshoot.energy import force_soil_temperature
from hydroshoot.hydraulic import soil_water_potential
from hydroshoot.params import Params


class HydroShootInputs(object):
    def __init__(self, path_project: Path, user_params: dict, scene: Scene,
                 is_write_result: bool = True, path_output_file: Path = None, **kwargs):
        self.path_project = path_project
        self.path_output_file = None
        self.path_output_dir = None
        self.scene = scene
        self.is_write_result = is_write_result

        if 'form_factors' in kwargs:
            self.set_form_factors(user_form_factors=kwargs['form_factors'])
        else:
            self.form_factors = None

        if 'leaf_nitrogen' in kwargs:
            self.leaf_nitrogen = self._parse_keys(user_input=kwargs['leaf_nitrogen'])
            self.is_nitrogen_calculated = True
        else:
            self.leaf_nitrogen = None
            self.is_nitrogen_calculated = False

        if 'leaf_ppfd' in kwargs:
            self.leaf_ppfd = self.set_ppfd(dicto=kwargs['leaf_ppfd'])
            self.is_ppfd_interception_calculated = True
        else:
            self.leaf_ppfd = None
            self.is_ppfd_interception_calculated = False

        self.psi_soil_forced = kwargs['psi_soil'] if 'psi_soil' in kwargs else None
        self.gdd_since_budbreak = kwargs['gdd_since_budbreak'] if 'gdd_since_budbreak' in kwargs else None
        self.sun2scene = kwargs['sun2scene'] if 'sun2scene' in kwargs else None
        self.soil_size = kwargs['soil_size'] if 'soil_size' in kwargs else None

        self.weather = None
        self.psi_pd = None

        self.params = Params(params_path=str(self.path_project / 'params.json'), user_params=user_params)
        self.set_weather()
        self.set_soil_predawn_water_potential()
        self.set_output_dir(user_path=path_output_file)

    def set_weather(self):
        df = read_csv(self.path_project / self.params.simulation.weather_file_name, sep=';', decimal='.', header=0)
        df.time = DatetimeIndex(df.time)
        df = df.set_index(df.time)
        if 'Ca' not in df.columns:
            df['Ca'] = [400.] * len(df)  # ppm [CO2]
        if 'Pa' not in df.columns:
            df['Pa'] = [101.3] * len(df)  # atmospheric pressure

        self.weather = df

    def set_soil_predawn_water_potential(self):
        if self.psi_soil_forced is None:
            assert (self.path_project / 'psi_soil.input').is_file(), "The 'psi_soil.input' file is missing."
            psi_pd = read_csv(self.path_project / 'psi_soil.input', sep=';', decimal='.').set_index('time')
            psi_pd.index = [datetime.strptime(str(s), "%Y-%m-%d") for s in psi_pd.index]
            self.psi_pd = psi_pd
        pass

    def set_output_dir(self, user_path: Path):
        if user_path is None:
            self.path_output_dir = self.path_project / 'output'
            self.path_output_file = self.path_output_dir / 'time_series.csv'
        else:
            self.path_output_file = user_path
            self.path_output_dir = self.path_output_file.parent

    def set_ppfd(self, dicto: dict) -> dict:
        res = {}
        for k1, v1 in dicto.items():
            res[k1] = {}
            for k, v in v1.items():
                if isinstance(v, dict):
                    res[k1][k] = self._parse_keys(user_input=v)
                else:
                    res[k1][k] = v
        return res

    def set_form_factors(self, user_form_factors: dict):
        self.form_factors = {k: self._parse_keys(user_input=v) for k, v in user_form_factors.items()}

    @staticmethod
    def _parse_keys(user_input: dict) -> dict:
        return {int(k): v for k, v in user_input.items()}


class HydroShootHourlyInputs(object):
    def __init__(self, psi_soil: float, sun2scene: Scene):
        self.date = None
        self.weather = None
        self.psi_soil = psi_soil
        self.sun2scene = sun2scene
        self.soil_temperature = None
        self.sky_temperature = None
        self.is_psi_soil_forced = self.psi_soil is not None

    def update(self, g: MTG, date_sim: datetime, hourly_weather: DataFrame, psi_pd: DataFrame, params: Params):
        self.date = date_sim
        self.weather = hourly_weather
        self.calc_psi_soil(g=g, psi_pd=psi_pd, params=params)
        self.sun2scene = visu(g, def_elmnt_color_dict=True, scene=Scene()) if self.sun2scene is not None else None
        self.soil_temperature = force_soil_temperature(self.weather)

        pass

    def calc_psi_soil(self, g: MTG, psi_pd: DataFrame, params: Params):
        if not self.is_psi_soil_forced:
            if self.date.hour == 0:
                try:
                    self.psi_soil = psi_pd.loc[self.date, :][0]
                except KeyError:
                    pass
            # Estimate soil water potential evolution due to transpiration
            else:
                self.psi_soil = soil_water_potential(
                    psi_soil_init=self.psi_soil,
                    water_withdrawal=g.node(g.node(g.root).vid_collar).Flux * params.simulation.conv_to_second,
                    soil_class=params.soil.soil_class,
                    soil_total_volume=params.soil.rhyzo_total_volume,
                    psi_min=params.hydraulic.psi_min)

        pass


def verify_inputs(g: MTG, inputs: HydroShootInputs):
    params = inputs.params

    # if inputs.psi_soil_forced is None:
    #     assert (inputs.path_project / 'psi_soil.input').is_file(), "The 'psi_soil.input' file is missing."

    assert params.irradiance.E_type in ('Rg_Watt/m2', 'RgPAR_Watt/m2' 'PPFD_umol/m2/s'), (
        "'irradiance_unit' must be one of the following 'Rg_Watt/m2', 'RgPAR_Watt/m2' or'PPFD_umol/m2/s'.")

    if 'length' in params.soil.soil_dimensions:
        assert params.soil.soil_dimensions['length'] <= params.planting.spacing_on_row, (
            "soil 'length' dimension must not exceed that of plant 'spacing_on_row'")
        assert params.soil.soil_dimensions['width'] <= params.planting.spacing_between_rows, (
            "soil 'width' dimension must not exceed that of plant 'spacing_between_rows'")
    elif 'radius' in params.soil.soil_dimensions:
        assert params.soil.soil_dimensions['radius'] <= 0.5 * min(
            [params.planting.spacing_on_row, params.planting.spacing_between_rows]), (
            "soil 'radius' dimension must not exceed the half of any of 'spacing_on_row' or 'spacing_between_rows'")

    if params.irradiance.E_type.split('_')[0] == 'PPFD':
        assert 'PPFD' in inputs.weather.columns, '"PPFD" column is missing in weather data'
    else:
        assert 'Rg' in inputs.weather.columns, '"Rg" column is missing in weather data'

    if not inputs.is_nitrogen_calculated:
        assert (params.simulation.date_beg - min(inputs.weather.index)).days >= 10, (
            'Meteorological data do not cover 10 days prior to simulation date.')
        if inputs.gdd_since_budbreak is None:
            assert min(inputs.weather.index) <= params.phenology.date_budbreak, (
                'The provided weather data do not cover the period from budbreak to the start date of simulation')

    if params.soil.rhyzo_solution:
        radius_max = 0.5 * min(params.soil.soil_dimensions['length'], params.soil.soil_dimensions['width'])
        assert max(params.soil.rhyzo_radii) <= radius_max, f'Maximum soil radius must not exceed {radius_max} cm'

    if inputs.leaf_nitrogen is not None:
        assert all([vid in inputs.leaf_nitrogen
                    for vid in get_leaves(g=g, leaf_lbl_prefix=inputs.params.mtg_api.leaf_lbl_prefix)])

    if inputs.leaf_ppfd is not None:
        for sim_datetime in inputs.params.simulation.date_range:
            d = sim_datetime.strftime("%Y%m%d%H%M%S")
            assert d in inputs.leaf_ppfd, f'Absorbed PPFD on {sim_datetime} is not provided'
            assert all([s in inputs.leaf_ppfd[d] for s in ('diffuse_to_total_irradiance_ratio', 'Ei', 'Eabs')])
            assert all([s in inputs.leaf_ppfd[d]['Ei'] for s in inputs.leaf_ppfd[d]['Eabs']])
            assert all([vid in inputs.leaf_ppfd[d]['Eabs']
                        for vid in get_leaves(g=g, leaf_lbl_prefix=inputs.params.mtg_api.leaf_lbl_prefix)])

    if inputs.form_factors is not None:
        assert all([s in inputs.form_factors for s in ('ff_sky', 'ff_leaves', 'ff_soil')]), (
            "At least one form factor key is missing or misspelt, "
            "('ff_sky', 'ff_leaves' and 'ff_soil' keys should be provided)")
        for k, v in inputs.form_factors.items():
            assert all([vid in v for vid in get_leaves(g=g, leaf_lbl_prefix=inputs.params.mtg_api.leaf_lbl_prefix)])

    if params.simulation.is_hydraulic_structure:
        gs_model = params.exchange.par_gs['model']
        assert gs_model != 'vpd', (
            f'The stomatal conductance model: "{gs_model}" does not require the hydraulic structure to be calculated')


def print_sim_infos(inputs: HydroShootInputs):
    params = inputs.params
    site_location = params.simulation.geo_location
    print("-" * 72)
    if inputs.gdd_since_budbreak is not None:
        print(f'GDD since budbreak: {inputs.gdd_since_budbreak} °Cd')
    print(f'Soil class: {params.soil.soil_class}')
    print(f'Site location: ({site_location[0]:.2f} N, {site_location[1]:.2f} E, {site_location[2]:.2f} m amsl)')
    print(f'Calculate Energy budget: {params.simulation.is_energy_budget}')
    print(f'Ignore shoot resistance: {params.simulation.is_negligible_shoot_resistance}')
    print(f'Calculate Hydraulic structure: {params.simulation.is_hydraulic_structure}')
    print(f'Add Rhyzoshpere cylinders: {params.soil.rhyzo_solution}')
    print(f'Use user form factors: {inputs.form_factors is not None}')
    print(f'Use user irradiance data: {inputs.leaf_ppfd is not None}')
    print(f'Use user nitrogen data: {inputs.leaf_nitrogen is not None}')
    print("-" * 72)
    pass
