from datetime import datetime

from pandas import date_range, DatetimeIndex, DataFrame


def calc_gdd_since_budbreak(weather: DataFrame, date_budbreak: datetime, date_beg_sim: datetime,
                            temperature_base: float) -> float:
    tdays = date_range(date_budbreak, date_beg_sim, freq='D')
    tmeteo = weather.loc[tdays, :].Tac.to_frame()
    tmeteo = tmeteo.set_index(DatetimeIndex(tmeteo.index).normalize())
    df_min = tmeteo.groupby(tmeteo.index).aggregate(min).Tac
    df_max = tmeteo.groupby(tmeteo.index).aggregate(max).Tac
    df_tt = 0.5 * (df_min + df_max) - temperature_base
    df = weather.resample('D').agg({'Tac': ['max', 'min']}).mean(axis=1)
    df.apply(lambda x: )
    return df_tt.cumsum()[-1]
