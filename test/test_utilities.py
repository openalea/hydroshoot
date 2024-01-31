from openalea.hydroshoot import utilities


def test_calc_effective_daily_temperature():
    assert (utilities.calc_effective_daily_temperature(temperature_air=10, temperature_base=10) == 0)
    assert (utilities.calc_effective_daily_temperature(temperature_air=9, temperature_base=10) == 0)
    assert (utilities.calc_effective_daily_temperature(temperature_air=11, temperature_base=10) == 1)
