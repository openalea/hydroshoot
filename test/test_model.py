""" A global test of hydroshoot model on potted grapevine, to secure refactoring"""
from pathlib import Path

from numpy.testing import assert_array_almost_equal

from hydroshoot import model
from test import non_regression_data


def test_potted_grapevine():
    results = model.run(
        g=non_regression_data.potted_syrah(),
        wd=Path(__file__).parent / 'data',
        write_result=False,
        psi_soil=-0.5,
        gdd_since_budbreak=1000.)
    ref = non_regression_data.reference_time_series_output()
    # do not compare date index
    assert_array_almost_equal(ref.iloc[0, 1:], results.reset_index(drop=True).iloc[0, :], decimal=0)
