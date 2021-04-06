import os
import sys
import unittest
import numpy as np

from octopuspy.info import Info

test_data_dir = 'test_data'
print(os.path.dirname(os.path.abspath(__file__)))
if 'win' in sys.platform:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\'
else:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '/'


class TestInfo(unittest.TestCase):
    ''' '''

    def test_filepath(self):
        '''Make sure the filepath is properly set'''

        info_data = Info(test_path)
        filepath = test_path + 'info'

        self.assertEqual(info_data._info_path, filepath)

    def test_num_ions(self):
        '''Make sure the number of ions is properly set'''

        info_data = Info(test_path)

        num_ions = 2
        self.assertEqual(info_data.num_ions, num_ions)

    def test_lattice_vector(self):
        '''Make sure the lattice vector is properly retrieved'''

        info_data = Info(test_path)

        zipped_vects = info_data.get_lattice_vectors()
        lattice_vector, reciprocal_lattice_vector = zip(*zipped_vects)
        lattice_vector = list(lattice_vector)

        lat_vect = ['    3.866975    0.000000    0.000000',
                    '    1.933488    3.348899    0.000000',
                    '    1.933488    1.116300    3.157372']

        self.assertEqual(lattice_vector, lat_vect)

    def test_reciprocal_lattice_vector(self):
        '''Make sure the reciprocal lattice vector is properly retrieved'''

        info_data = Info(test_path)

        zipped_vects = info_data.get_lattice_vectors()
        lattice_vector, reciprocal_lattice_vector = zip(*zipped_vects)
        reciprocal_lattice_vector = list(reciprocal_lattice_vector)

        recp_lat_vect = ['    1.624832   -0.938097   -0.663335',
                         '    0.000000    1.876195   -0.663335',
                         '    0.000000    0.000000    1.990005']

        self.assertEqual(reciprocal_lattice_vector, recp_lat_vect)

if __name__ == '__main__':
    unittest.main()
