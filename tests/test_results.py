import os
import sys
import unittest
import numpy as np

from octopuspy.results import Results

test_data_dir = 'test_data'
if 'win' in sys.platform:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\'
else:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '/'


class TestResults(unittest.TestCase):
    ''' '''

    def test_filepath(self):
        '''Make sure the filepath is properly set'''

        num_kpoints = 4
        results_data = Results(test_path, num_kpoints)
        filepath = test_path + 'results.out'

        self.assertEqual(results_data._results_path, filepath)

    def test_weights(self):
        '''Make sure the weights are properly extracted from the results.out file'''

        num_kpoints = 4
        results_data = Results(test_path, num_kpoints)

        weights = results_data.weights

        test_weights = np.array(['0.004630', '0.004630', '0.004630', '0.004630'])
        test_weights = test_weights.astype('float64')

        self.assertTrue(np.allclose(weights, test_weights))


if __name__ == '__main__':
    unittest.main()
