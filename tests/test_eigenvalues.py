import os
import sys
import unittest
import numpy as np

from octopuspy.eigenvalues import Eigenvalues

test_data_dir = 'test_data'
if 'win' in sys.platform:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\'
    partials_path = test_path + '\\partial_occupancies\\'
else:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '/'
    partials_path = test_path + '/partial_occupancies/'


class TestEigenvalues(unittest.TestCase):
    '''Many tests are redundant with test_bandstructure.py however this covers partial occupancies '''

    def test_get_occupancies_no_partials(self):
        '''Tests the setting of the occuancies with no partial occupancies'''

        eigenv = Eigenvalues(test_path, num_kpoints=4, num_bands=5)

        occupancies =  np.array([[2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0]])
        eignv_occ = eigenv.get_occupancies()
        self.assertTrue(np.allclose(eignv_occ, occupancies))

    def test_get_occupancies_with_partials(self):
        '''Tests the setting of the occuancies with partial occupancies'''

        eigenv = Eigenvalues(partials_path, num_kpoints=4, num_bands=5)

        occupancies =  np.array([[2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0],
                                 [2.0, 2.0, 2.0, 2.0, 0.0]])
        eignv_occ = eigenv.get_occupancies()
        self.assertTrue(np.allclose(eignv_occ, occupancies))


if __name__ == '__main__':
    unittest.main()