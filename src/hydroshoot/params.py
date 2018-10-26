#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:10:00 2018

@author: rami
"""

from json import load
from jsonschema import validate, validators

json_file_path = r'C:\Users\ralbasha\Documents\dvlp\hydroshoot\src\hydroshoot_data\params.json'
json_schm_path = r'C:\Users\ralbasha\Documents\dvlp\hydroshoot\src\hydroshoot_data\params_schema.json'
json_file = load(open(json_file_path, mode='r', encoding="utf-8"))
json_schm = load(open(json_schm_path, mode='r', encoding="utf-8"))
validate(json_file, json_schm)