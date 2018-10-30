#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:10:00 2018

@author: rami
"""

import os
from json import load
from jsonschema import validate, validators


class Params():

    def __init__(self, params_path):
        self.params_path = params_path
        self.params_schema = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '../hydroshoot/params_schema.json'))

        user_params = self._get_user_params()

        self.simulation = Simulation(user_params['simulation'])
        self.phenology = Phenology(user_params['phenology'])
        self.mtg_api = MTG_API(user_params['mtg_api'])
        self.numerical_resolution = NumericalResolution(user_params['numerical_resolution'])


    def _get_user_params(self):
        """
        Get parameters values defined by the user.

        :Params:
            - **params_path**: (string) absolute path to parameters json file

        :Return:
            - (dict) dictionary of parameters identifiers and values.
        """

        json_file = load(open(self.params_path, mode='r', encoding="utf-8"))
        json_schm = load(open(json_schm_path, mode='r', encoding="utf-8"))
        validate(json_file, json_schm)

        return json_file

class Simulation():

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


class Phenology:

    def __init__(self, phenology_dict):
        self.emdate = phenology_dict['emdate']
        self.t_base = phenology_dict['t_base']


class MTG_API():

    def __init__(self, mtg_dict):
        self.collar_label = mtg_dict['collar_label']
        self.leaf_lbl_prefix = mtg_dict['leaf_lbl_prefix']
        self.stem_lbl_prefix = mtg_dict['stem_lbl_prefix']


class NumericalResolution():

    def __init__(self, numerical_resolution_dict):
        self.max_iter = numerical_resolution_dict['max_iter']
        self.psi_step = numerical_resolution_dict['psi_step']
        self.psi_error_threshold = numerical_resolution_dict['psi_error_threshold']
        self.t_step = numerical_resolution_dict['t_step']
        self.t_error_crit = numerical_resolution_dict['t_error_crit']




    #     "irradiance": {
    #         "E_type": "Rg_Watt/m2",
    #         "E_type2": "Ei",
    #         "opt_prop": {
    #             "SW": {
    #                 "leaf": [
    #                     0.06,
    #                     0.07
    #                 ],
    #                 "stem": [
    #                     0.13
    #                 ],
    #                 "other": [
    #                     0.06,
    #                     0.07
    #                 ]
    #             },
    #             "LW": {
    #                 "leaf": [
    #                     0.04,
    #                     0.07
    #                 ],
    #                 "stem": [
    #                     0.13
    #                 ],
    #                 "other": [
    #                     0.06,
    #                     0.07
    #                 ]
    #             }
    #         },
    #         "scene_rotation": 0.0,
    #         "turtle_format": "soc",
    #         "turtle_sectors": "46",
    #         "icosphere_level": "None"
    #     },
    #     "energy": {
    #         "solo": true,
    #         "limit": -0.000000001,
    #         "t_cloud": 2.0,
    #         "t_sky": -20.0
    #     },
    #     "hydraulic": {
    #         "MassConv": 18.01528,
    #         "psi_min": -3.0,
    #         "Kx_dict": {
    #             "a": 1.6,
    #             "b": 2.0,
    #             "min_kmax": 0.01
    #         },
    #         "par_K_vul": {
    #             "model": "misson",
    #             "fifty_cent": -0.51,
    #             "sig_slope": 1.0
    #         }
    #     },
    #     "exchange": {
    #         "rbt": 0.6667,
    #         "ca": 360.0,
    #         "Na_dict": {
    #             "aN": -0.0008,
    #             "bN": 3.3,
    #             "aM": 6.471,
    #             "bM": 56.635
    #         },
    #         "par_gs": {
    #             "model": "vpd",
    #             "g0": 0.0,
    #             "m0": 5.278,
    #             "psi0": -1.0,
    #             "D0": 30.0,
    #             "n": 4.0
    #         },
    #         "par_photo": {
    #             "Vcm25": 89.0,
    #             "Jm25": 143.0,
    #             "cRd": 0.008,
    #             "TPU25": 10.0,
    #             "Kc25": 404.9,
    #             "Ko25": 278.4,
    #             "Tx25": 42.75,
    #             "alpha": [
    #                 0.2,
    #                 0.2,
    #                 0.19,
    #                 0.19,
    #                 0.14,
    #                 0.12
    #             ],
    #             "alpha_T_limit": [
    #                 15,
    #                 20,
    #                 25,
    #                 30,
    #                 34,
    #                 50
    #             ],
    #             "a1": 0.98,
    #             "a2": 0.98,
    #             "a3": 0.98,
    #             "ds": 0.635,
    #             "dHd": 200.0,
    #             "RespT_Kc": {
    #                 "model": "Arrhenius",
    #                 "c": 38.05,
    #                 "deltaHa": 79.43
    #             },
    #             "RespT_Ko": {
    #                 "model": "Arrhenius",
    #                 "c": 20.30,
    #                 "deltaHa": 36.38
    #             },
    #             "RespT_Vcm": {
    #                 "model": "Arrhenius",
    #                 "c": 26.35,
    #                 "deltaHa": 65.33
    #             },
    #             "RespT_Jm": {
    #                 "model": "Arrhenius",
    #                 "c": 17.57,
    #                 "deltaHa": 43.54
    #             },
    #             "RespT_TPU": {
    #                 "model": "Arrhenius",
    #                 "c": 21.46,
    #                 "deltaHa": 53.1
    #             },
    #             "RespT_Rd": {
    #                 "model": "Arrhenius",
    #                 "c": 18.72,
    #                 "deltaHa": 46.39
    #             },
    #             "RespT_Tx": {
    #                 "model": "Arrhenius",
    #                 "c": 19.02,
    #                 "deltaHa": 37.83
    #             }
    #         },
    #         "par_photo_N": {
    #             "Vcm25_N": [
    #                 34.02,
    #                 -3.13
    #             ],
    #             "Jm25_N": [
    #                 78.27,
    #                 -17.3
    #             ],
    #             "Rd_N": [
    #                 0.42,
    #                 -0.01
    #             ],
    #             "TPU25_N": [
    #                 6.24,
    #                 -1.92
    #             ]
    #         }
    #     },
    #     "soil": {
    #         "soil_class": "Sand",
    #         "soil_dimensions": [
    #             3.6,
    #             1.0,
    #             1.2
    #         ],
    #         "rhyzo_solution": true,
    #         "rhyzo_number": 3,
    #         "rhyzo_radii": [
    #             0.167,
    #             0.333,
    #             0.5
    #         ],
    #         "rhyzo_coeff": 0.5,
    #         "roots": [
    #             0.013,
    #             0.0001
    #         ]
    #     }
    # }
    #
