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
