#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:10:00 2018

@author: rami
"""

import os
from copy import deepcopy
from datetime import datetime
from json import load
from math import pi

from jsonschema import validate
from pandas import date_range


class Params:

    def __init__(self, params_path: str = None, user_params: dict = None):
        if user_params is None:
            user_params = self._get_user_params(path_params=params_path)

        self._validate_user_params(
            user_params=user_params,
            path_params_schema=os.path.abspath(
                os.path.join(os.path.dirname(__file__), '../hydroshoot/params_schema.json')))
        self.simulation = Simulation(user_params['simulation'])
        self.phenology = Phenology(user_params['phenology'])
        self.mtg_api = MtgAPI(user_params['mtg_api'])
        self.numerical_resolution = NumericalResolution(user_params['numerical_resolution'])
        self.irradiance = Irradiance(user_params['irradiance'])
        self.energy = Energy(user_params['energy'])
        self.hydraulic = Hydraulic(user_params['hydraulic'])
        self.exchange = Exchange(user_params['exchange'])
        self.planting = Planting(user_params['planting'])
        self.soil = Soil(soil_dict=user_params['soil'], planting_dict=user_params['planting'])

        self._set_pattern()

    @staticmethod
    def _get_user_params(path_params: str) -> str:
        """Gets parameter values defined by the user.

        Args:
            path_params: absolute path to parameters json file

        Returns:
            - model parameters
        """

        with open(path_params, mode='r') as f:
            params = load(f)

        return params

    @staticmethod
    def _validate_user_params(user_params: dict, path_params_schema: str) -> None:
        """Validate parameter values defined by the user.

        Args:
            user_params: model parameters defined by the user

        """

        with open(path_params_schema, mode='r') as f:
            params_scehma = load(f)

        validate(user_params, params_scehma)
        pass

    def _set_pattern(self):
        y_max = self.planting.spacing_between_rows / self.simulation.conv_to_meter
        x_max = self.planting.spacing_on_row / self.simulation.conv_to_meter
        self.irradiance.pattern = ((-x_max / 2.0, -y_max / 2.0), (x_max / 2.0, y_max / 2.0))


class Simulation:

    def __init__(self, simulation_dict):
        self.date_beg = datetime.strptime(simulation_dict['sdate'], "%Y-%m-%d %H:%M:%S")
        self.date_end = datetime.strptime(simulation_dict['edate'], "%Y-%m-%d %H:%M:%S")
        self._latitude = simulation_dict['latitude']
        self._longitude = simulation_dict['longitude']
        self._elevation = simulation_dict['elevation']
        self.tzone = simulation_dict['tzone']
        self.output_index = simulation_dict['output_index']
        self.unit_scene_length = simulation_dict['unit_scene_length']
        self.is_hydraulic_structure = simulation_dict['hydraulic_structure']
        self.is_negligible_shoot_resistance = simulation_dict['negligible_shoot_resistance']
        self.is_energy_budget = simulation_dict['energy_budget']
        self.weather_file_name = simulation_dict['meteo']

        self.date_range = date_range(start=self.date_beg, end=self.date_end, freq='H')
        self.conv_to_second = {'D': 86.4e3, 'H': 3600., 'T': 60., 'S': 1.}[self.date_range.freqstr]
        self.conv_to_meter = {'mm': 1.e-3, 'cm': 1.e-2, 'm': 1.}[self.unit_scene_length]
        self.geo_location = (self._latitude, self._longitude, self._elevation)


class Phenology:

    def __init__(self, phenology_dict):
        self.date_budbreak = datetime.strptime(phenology_dict['emdate'], "%Y-%m-%d %H:%M:%S")
        self.t_base = phenology_dict['t_base']


class MtgAPI:

    def __init__(self, mtg_dict):
        self.collar_label = mtg_dict['collar_label']
        self.leaf_lbl_prefix = mtg_dict['leaf_lbl_prefix']
        self.stem_lbl_prefix = tuple(mtg_dict['stem_lbl_prefix'])


class NumericalResolution:

    def __init__(self, numerical_resolution_dict):
        self.max_iter = numerical_resolution_dict['max_iter']
        self.psi_step = numerical_resolution_dict['psi_step']
        self.psi_error_threshold = numerical_resolution_dict['psi_error_threshold']
        self.t_step = numerical_resolution_dict['t_step']
        self.t_error_threshold = numerical_resolution_dict['t_error_threshold']


class Planting:
    def __init__(self, planting_dict):
        self.spacing_between_rows = planting_dict['spacing_between_rows']
        self.spacing_on_row = planting_dict['spacing_on_row']
        self.scene_rotation = planting_dict['row_angle_with_south']


class Irradiance:

    def __init__(self, irradiance_dict):
        self.E_type = irradiance_dict['E_type']
        self.E_type2 = irradiance_dict['E_type2']
        self.opt_prop = _list2tuple(irradiance_dict['opt_prop'])
        self.turtle_format = irradiance_dict['turtle_format']
        self.turtle_sectors = irradiance_dict['turtle_sectors']
        self.icosphere_level = irradiance_dict['icosphere_level']

        self.pattern = None


class Energy:

    def __init__(self, energy_dict):
        self.solo = energy_dict['solo']
        self.t_cloud = energy_dict['t_cloud']
        self.t_sky = energy_dict['t_sky']


class Hydraulic:

    def __init__(self, hydraulic_dict):
        self.psi_min = hydraulic_dict['psi_min']
        self.Kx_dict = hydraulic_dict['Kx_dict']
        self.par_K_vul = hydraulic_dict['par_K_vul']


class Exchange:

    def __init__(self, exchange_dict):
        self.rbt = exchange_dict['rbt']
        self.Na_dict = exchange_dict['Na_dict']
        self.par_gs = exchange_dict['par_gs']
        self.par_photo = {
            "photo_inhibition": {
                "dhd_inhib_beg": 195,
                "dHd_inhib_max": 180,
                "psi_inhib_beg": -0.75,
                "psi_inhib_max": -2,
                "temp_inhib_beg": 32,
                "temp_inhib_max": 33}}
        self.par_photo.update(exchange_dict['par_photo'])
        self.par_photo_N = exchange_dict['par_photo_N']


class Soil:

    def __init__(self, soil_dict, planting_dict):
        self.soil_class = soil_dict['soil_class']
        self.soil_dimensions = soil_dict['soil_dimensions']
        self.rhyzo_coeff = soil_dict['rhyzo_coeff']

        if all([s in soil_dict for s in ('rhyzo_coeff', 'avg_root_radius', 'root_length')]):
            self.rhyzo_solution = True
            self.avg_root_radius = soil_dict['avg_root_radius']
            self.root_length = soil_dict['avg_root_radius']
            self.rhyzo_coeff = soil_dict['rhyzo_coeff']
        else:
            self.rhyzo_solution = False
            self.avg_root_radius = None
            self.root_length = None
            self.rhyzo_coeff = 1

        soil_dimensions = self.get_soil_dimensions(soil_data=soil_dict['soil_dimensions'], planting_data=planting_dict)
        self.soil_volume = self.calc_soil_volume(soil_dimensions=soil_dimensions)
        self.rhyzo_volume = self.rhyzo_coeff * self.soil_volume

    @staticmethod
    def calc_soil_volume(soil_dimensions: dict) -> float:
        if 'radius' in soil_dimensions:
            res = pi * (soil_dimensions['radius'] ** 2) * soil_dimensions['depth']
        else:
            res = soil_dimensions['length'] * soil_dimensions['width'] * soil_dimensions['depth']
        return res

    @staticmethod
    def get_soil_dimensions(soil_data: dict, planting_data: dict) -> dict:
        if 'radius' not in soil_data and not all([s in soil_data for s in ('width', 'length')]):
            res = {
                'width': planting_data['spacing_between_rows'],
                'length': planting_data['spacing_on_row'],
                'depth': soil_data['depth']}
        else:
            res = soil_data
        return res


def _list2tuple(dct):
    _dct = deepcopy(dct)
    for k_wave, v_wave in dct.items():
        for k, v in v_wave.items():
            _dct[k_wave][k] = tuple(v)
    return _dct
