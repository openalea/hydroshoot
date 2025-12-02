""" A global test of hydroshoot model on potted grapevine, to secure refactoring"""
from pathlib import Path

from numpy.testing import assert_array_almost_equal

from openalea.hydroshoot import model
import non_regression_data


def test_potted_grapevine():
    path_data = Path(__file__).parent / 'data'
    results = model.run(
        g=non_regression_data.potted_syrah(),
        wd=path_data,
        path_weather=path_data / 'meteo.csv',
        write_result=False,
        psi_soil=-0.2,
        gdd_since_budbreak=100.)
    ref = non_regression_data.reference_time_series_output()
    # do not compare date index
    assert_array_almost_equal(ref.iloc[0, :], results.reset_index(drop=True).iloc[0, :], decimal=0)
