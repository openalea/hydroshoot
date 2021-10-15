#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:10:00 2018

@author: rami
"""

import os
from copy import deepcopy
from json import load
from jsonschema import validate


class Params:

    def __init__(self, params_path):
        self._params_path = params_path
        self._params_schema = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '../hydroshoot/params_schema.json'))

        user_params = self._get_user_params()

        self.simulation = Simulation(user_params['simulation'])
        self.phenology = Phenology(user_params['phenology'])
        self.mtg_api = MtgAPI(user_params['mtg_api'])
        self.numerical_resolution = NumericalResolution(user_params['numerical_resolution'])
        self.irradiance = Irradiance(user_params['irradiance'])
        self.energy = Energy(user_params['energy'])
        self.hydraulic = Hydraulic(user_params['hydraulic'])
        self.exchange = Exchange(user_params['exchange'])
        self.soil = Soil(user_params['soil'])

    def _get_user_params(self):
        """
        Get parameters values defined by the user.

        :Params:
            - **params_path**: (string) absolute path to parameters json file

        :Return:
            - (dict) dictionary of parameters identifiers and values.
        """

        with open(self._params_path, mode='r') as f:
            json_file = load(f, encoding="utf-8")

        with open(self._params_schema, mode='r') as f:
            json_schm = load(f, encoding="utf-8")

        validate(json_file, json_schm)

        return json_file


class Simulation:

    def __init__(self, simulation_dict):
        self.sdate = simulation_dict['sdate']
        self.edate = simulation_dict['edate']
        self.latitude = simulation_dict['latitude']
        self.longitude = simulation_dict['longitude']
        self.elevation = simulation_dict['elevation']
        self.tzone = simulation_dict['tzone']
        self.output_index = simulation_dict['output_index']
        self.unit_scene_length = simulation_dict['unit_scene_length']
        self.simplified_form_factors = simulation_dict['simplified_form_factors']
        self.hydraulic_structure = simulation_dict['hydraulic_structure']
        self.negligible_shoot_resistance = simulation_dict['negligible_shoot_resistance']
        self.energy_budget = simulation_dict['energy_budget']
        self.soil_water_deficit = simulation_dict['soil_water_deficit']
        self.meteo = simulation_dict['meteo']


class Phenology:

    def __init__(self, phenology_dict):
        self.emdate = phenology_dict['emdate']
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
        self.t_error_crit = numerical_resolution_dict['t_error_crit']


class Irradiance:

    def __init__(self, irradiance_dict):
        self.E_type = irradiance_dict['E_type']
        self.E_type2 = irradiance_dict['E_type2']
        self.opt_prop = _list2tuple(irradiance_dict['opt_prop'])
        self.scene_rotation = irradiance_dict['scene_rotation']
        self.turtle_format = irradiance_dict['turtle_format']
        self.turtle_sectors = irradiance_dict['turtle_sectors']
        self.icosphere_level = irradiance_dict['icosphere_level']


class Energy:

    def __init__(self, energy_dict):
        self.solo = energy_dict['solo']
        self.limit = energy_dict['limit']
        self.t_cloud = energy_dict['t_cloud']
        self.t_sky = energy_dict['t_sky']


class Hydraulic:

    def __init__(self, hydraulic_dict):
        self.MassConv = hydraulic_dict['MassConv']
        self.psi_min = hydraulic_dict['psi_min']
        self.Kx_dict = hydraulic_dict['Kx_dict']
        self.par_K_vul = hydraulic_dict['par_K_vul']


class Exchange:

    def __init__(self, exchange_dict):
        self.rbt = exchange_dict['rbt']
        self.ca = exchange_dict['ca']
        self.Na_dict = exchange_dict['Na_dict']
        self.par_gs = exchange_dict['par_gs']
        self.par_photo = exchange_dict['par_photo']
        self.par_photo_N = exchange_dict['par_photo_N']


class Soil:

    def __init__(self, soil_dict):
        self.soil_class = soil_dict['soil_class']
        self.soil_dimensions = soil_dict['soil_dimensions']
        self.rhyzo_solution = soil_dict['rhyzo_solution']
        self.rhyzo_radii = soil_dict['rhyzo_radii']
        self.rhyzo_coeff = soil_dict['rhyzo_coeff']
        self.roots = soil_dict['roots']


def _list2tuple(dct):
    _dct = deepcopy(dct)
    for k_wave, v_wave in dct.items():
        for k, v in v_wave.items():
            _dct[k_wave][k] = tuple(v)
    return _dct
