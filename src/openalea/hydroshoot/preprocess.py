from datetime import datetime

from pandas import date_range, DataFrame

from openalea.hydroshoot.utilities import calc_effective_daily_temperature


def calc_gdd_since_budbreak(weather: DataFrame, date_budbreak: datetime, date_beg_sim: datetime,
                            temperature_base: float) -> float:
    time_range = date_range(date_budbreak, date_beg_sim, freq='H')
    df = weather.loc[time_range, :].resample('D').agg({'Tac': ['max', 'min']}).mean(axis=1)
    return df.apply(lambda x:
                    calc_effective_daily_temperature(temperature_air=x, temperature_base=temperature_base)).cumsum()[-1]
