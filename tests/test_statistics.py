"""Small numerical checks using invented values, never clinical records."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from reanalyze import bh, quantile_normalize, exact_mean_permutation

class StatisticsTests(unittest.TestCase):
    def test_bh_unsorted_and_missing(self):
        actual=bh([.04,.01,np.nan,.03])
        np.testing.assert_allclose(actual[[0,1,3]],[.04,.03,.04])
        self.assertTrue(np.isnan(actual[2]))

    def test_ties_average_reference_positions(self):
        # Reference quantiles are [1, 2, 4]. First column's tied zeros
        # occupy positions 1 and 2, so the corrected value is 1.5.
        data=pd.DataFrame({'a':[0.,0.,2.],'b':[2.,4.,6.]})
        np.testing.assert_allclose(quantile_normalize(data,'average_ties')['a'],[1.5,1.5,4.])
        np.testing.assert_allclose(quantile_normalize(data)['a'],[1.,1.,4.])

    def test_exact_two_sided_permutation(self):
        # Of the six 2-vs-2 allocations, the two opposite extreme
        # partitions are at least as extreme as the observed one.
        data=pd.DataFrame([[0.,1.,2.,3.]])
        np.testing.assert_allclose(exact_mean_permutation(data,np.array([True,True,False,False])),[2/6])

if __name__=='__main__':unittest.main()
