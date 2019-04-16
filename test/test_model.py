""" A global test of hydroshoot model on potted grapevine, to secure refactoring"""
import os
from hydroshoot.data_access import get_data_dir
from hydroshoot.data import potted_syrah, reference_time_series_output
from hydroshoot import model
from numpy.testing import assert_array_almost_equal


def test_potted_grapevine():
    g = potted_syrah()
    dir = os.path.join(get_data_dir(), '')
    results = model.run(g, dir, write_result=False,psi_soil=-0.5,
              gdd_since_budbreak=1000.)
    ref = reference_time_series_output()
    #do not compare date index
    assert_array_almost_equal(ref.iloc[0, 1:], results.reset_index(drop=True).iloc[0, :])