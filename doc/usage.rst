=====
Usage
=====

Set the required files up
=========================

1. Model parameters (fixed values)

Almost all model parameters are setup in the **params.json** file.
*params.json* is a *JavaScript Object Notation (JSON)* file that contains the following categories:

*  "simulation"
*  "phenology"
*  "mtg_api"
*  "numerical_resolution"
*  "irradiance"
*  "energy"
*  "hydraulic"
*  "exchange"
*  "soil"

Each of the above catgories referes to its hyponyme related processes. An exaustive and thoroughly detailed desctiption
of all the parameters inside each categries, their units and their expected values are given in the `params_schema.json`
file (src/hydroshoot/). The reader is encouraged to refer to this file to know how exactly to fill up the "params.json"
file.

2. Meteorological data

Meteorological data must be provided using a *comma-separated values* (csv) file whose name must be given in the
`params.json` file ("meteo" parameters)
This file must contains the following columns:
    * `time`: a `datetime` string having the format YYYY-mm-DD HH:MM:SS
    * `Tac`: actual air temperature as measured by the meteorological station :math:`[^\circ C]`
    * `hs`: relative humidity (%)
    * `u`: wind speed :math:`[m \ s^{-1}]`
    * `Rg`: solar radiation (shortwave irradiance) given even in :math:`[W_{global} \ m_{ground}^{-2}]`,
      :math:`[W_{PAR} \ m_{ground}^{-2}]` or :math:`[\mu mol \ m_{ground}^{-2} \ s^{-1}]`. The user must provide the
      unit of `Rg` in the `params.json` file ("E_type" parameter).
