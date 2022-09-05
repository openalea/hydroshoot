from datetime import datetime

from pandas import DataFrame, read_csv, DatetimeIndex

from hydroshoot.params import Params


def init_weather(path_project: str, params: Params) -> DataFrame:
    meteo_tab = read_csv(f'{path_project}{params.simulation.meteo}', sep=';', decimal='.', header=0)
    meteo_tab.time = DatetimeIndex(meteo_tab.time)
    meteo_tab = meteo_tab.set_index(meteo_tab.time)

    #   Adding missing data
    if 'Ca' not in meteo_tab.columns:
        meteo_tab['Ca'] = [400.] * len(meteo_tab)  # ppm [CO2]
    if 'Pa' not in meteo_tab.columns:
        meteo_tab['Pa'] = [101.3] * len(meteo_tab)  # atmospheric pressure

    return meteo_tab


def init_soil_predawn_water_potential(path_project: str, params: Params, **kwargs) -> DataFrame:
    if 'psi_soil' in kwargs:
        psi_pd = DataFrame({'psi': [kwargs['psi_soil']] * len(params.simulation.date_range)},
                           index=params.simulation.date_range)
    else:
        psi_pd = read_csv(f'{path_project}psi_soil.input', sep=';', decimal='.').set_index('time')
        psi_pd.index = [datetime.strptime(str(s), "%Y-%m-%d") for s in psi_pd.index]

    return psi_pd
