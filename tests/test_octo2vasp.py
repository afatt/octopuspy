import os
import re
import sys
import unittest
import numpy as np
from glob import glob
from io import StringIO

from octopuspy.octo2vasp import Octo2Vasp
from octopuspy.bandstructure import Bandstructure

test_data_dir = 'test_data'
if 'win' in sys.platform:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\'
else:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '/'

class TestOcto2Vasp(unittest.TestCase):
    ''' '''

    def test_gen_procar_pos_kpoints(self):
        '''
        Test Case 1: positive kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('./**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '.\\tests\\test_data\\pos_kpoints\\bandstructure'
        lin_path = './tests/test_data/pos_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', energy_scale=1.0, valence_band_index=4)

        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')

        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                # print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)

    def test_gen_procar_neg_kpoints(self):
        '''
        Test Case 2: negative kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('./**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '.\\tests\\test_data\\neg_kpoints\\bandstructure'
        lin_path = './tests/test_data/neg_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', energy_scale=1.0, valence_band_index=4)

        o2v.filepath = './test_data'
        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')
        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                # print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)

    def test_gen_procar_mixed_sign_kpoints(self):
        '''
        Test Case 3: mixed sign kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('./**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '.\\tests\\test_data\\mixed_sign_kpoints\\bandstructure'
        lin_path = './tests/test_data/mixed_sign_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', energy_scale=1.0, valence_band_index=4)

        o2v.filepath = './test_data'
        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')
        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)


class TestBandstructure(unittest.TestCase):
    ''' '''

    def test_name(self):
        '''Make sure the name attibute is properly set'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        name = 'Si_03082021'
        self.assertEqual(bs._name, name)

    def test_filepath(self):
        '''Make sure the filepath attibute is properly set'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        filepath = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\bandstructure'
        self.assertEqual(bs._bandstructure_path, filepath)

    def test_bandstructure(self):
        '''Make sure the bandstructure attibute is properly set'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        test_bs = np.array([[0.00000000, 0.00000000, 0.00000000, 0.00000000, -7.74445402, 4.07037476, 4.07037593, 4.07038233, 6.61983340],
                            [0.01373858, 0.04166667, 0.00000000, 0.04166667, -7.71460976, 3.86361199, 3.94921615, 3.94921884, 6.54685076],
                            [0.02747716, 0.08333333, 0.00000000, 0.08333333, -7.62519023, 3.35857616, 3.65323831, 3.65323944, 6.34969314],
                            [0.04121574, 0.12500000, 0.00000000, 0.12500000, -7.47653932, 2.71474520, 3.28609405, 3.28609571, 6.07653663]])

        self.assertTrue(np.allclose(bs._bandstructure, test_bs))

    def test_efermi_path(self):
        '''Tests the setting of the efermi filepath'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        filepath = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\total-dos-efermi.dat'

        self.assertEqual(bs._efermi_path, filepath)

    def test_kpoints(self):
        '''Tests the setting of the kpoints'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        test_kx = np.array([0.00000000, 0.04166667, 0.08333333, 0.12500000])
        test_ky = np.array([0.00000000, 0.00000000, 0.00000000, 0.00000000])
        test_kz = np.array([0.00000000, 0.04166667, 0.08333333, 0.12500000])

        kx, ky, kz = zip(*bs.kpoints)

        self.assertTrue(np.allclose(kx, test_kx))
        self.assertTrue(np.allclose(ky, test_ky))
        self.assertTrue(np.allclose(kz, test_kz))

    def test_num_bands(self):
        '''Tests the setting of the number of bands attribute'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        num_bands = 5

        self.assertEqual(bs.num_bands, num_bands)

    def test_num_kpoints(self):
        '''Tests the setting of the number of k-points attribute'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        num_kpoints = 4

        self.assertEqual(bs.num_kpoints, num_kpoints)

    def test_efermi(self):
        '''Tests the loading of the fermi energy'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        efermi = 4.070390
        self.assertEqual(bs.efermi, efermi)

    def test_energies(self):
        '''Tests the loading of the energies from the bandstructure'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)

        energies = np.array([[-7.74445402, 4.07037476, 4.07037593, 4.07038233, 6.61983340],
                             [-7.71460976, 3.86361199, 3.94921615, 3.94921884, 6.54685076],
                             [-7.62519023, 3.35857616, 3.65323831, 3.65323944, 6.34969314],
                             [-7.47653932, 2.71474520, 3.28609405, 3.28609571, 6.07653663]])

        energies = energies - 4.070390

        bs_e, bs_o = bs.get_eigenvalues()
        self.assertTrue(np.allclose(bs_e, energies))

    def test_occupancies(self):
        '''Tests the setting of the occupancies'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)
        occupancies =  np.array([[2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0]])
        bs_e, bs_o = bs.get_eigenvalues()
        self.assertTrue(np.allclose(bs_o, occupancies))

    def test_occupied_bands(self):
        '''Tests the setting of the occupied bands from the _split_bands method'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)

        occupied_bands = np.array([[-7.74445402, -7.71460976, -7.62519023, -7.47653932],
                                   [4.07037476, 3.86361199, 3.35857616, 2.71474520],
                                   [4.07037593, 3.94921615, 3.65323831, 3.28609405],
                                   [4.07038233, 3.94921884, 3.65323944, 3.28609571]])

        occupied_bands = occupied_bands - 4.070390
        bs_e, bs_o = bs.get_eigenvalues()
        oc, unoc = bs._split_bands()
        self.assertTrue(np.allclose(oc, occupied_bands))

    def test_unoccupied_bands(self):
        '''Tests the setting of the unoccupied bands from the _split_bands method'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)

        unoccupied_bands = np.array([[6.61983340, 6.54685076, 6.34969314, 6.07653663]])
        unoccupied_bands = unoccupied_bands - 4.070390
        bs_e, bs_o = bs.get_eigenvalues()
        oc, unoc = bs._split_bands()
        self.assertTrue(np.allclose(unoc, unoccupied_bands))

    def test_conduction_band(self):
        '''Tests the setting of the conduction band'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)

        # there is only one unoccupied band which also makes it the valence band
        unoccupied_bands = np.array([[6.61983340, 6.54685076, 6.34969314, 6.07653663]])
        band_min = 6.07653663

        cb, cb_min = bs._get_conduction_band(unoccupied_bands)
        self.assertTrue(np.allclose(cb, unoccupied_bands))
        self.assertEqual(cb_min, band_min)

    def test_valence_band(self):
        '''Tests the setting of the valence band'''

        bs = Bandstructure('Si_03082021', test_path, 1.0, 4)

        occupied_bands = np.array([[-7.74445402, -7.71460976, -7.62519023, -7.47653932],
                                   [4.07037476, 3.86361199, 3.35857616, 2.71474520],
                                   [4.07037593, 3.94921615, 3.65323831, 3.28609405],
                                   [4.07038233, 3.94921884, 3.65323944, 3.28609571]])

        val_band = occupied_bands[-1,:]
        band_max = 4.07038233

        vb, vb_max = bs._get_valence_band(occupied_bands)
        self.assertTrue(np.allclose(vb, val_band))
        self.assertEqual(vb_max, band_max)

if __name__ == '__main__':
    unittest.main()
