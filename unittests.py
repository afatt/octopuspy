# procar testing test negative and positive kpoint values for x, y, and z
# test against the regular expression

import re
import unittest
import numpy as np

import bandstructure

# class TestGenProcar(unittest.TestCase):
#     ''' '''


class TestBandstructure(unittest.TestCase):
    ''' '''

    def test_name(self):
        '''Make sure the name attibute is properly set'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        name = 'Si_03082021'
        self.assertEqual(bs._name, name)

    def test_filepath(self):
        '''Make sure the filepath attibute is properly set'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        filepath = './ut_data/bandstructure'
        self.assertEqual(bs._bandstructure_path, filepath)

    def test_bandstructure(self):
        '''Make sure the bandstructure attibute is properly set'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        test_bs = np.array([[0.00000000, 0.00000000, 0.00000000, 0.00000000, -7.74445402, 4.07037476, 4.07037593, 4.07038233, 6.61983340],
                            [0.01373858, 0.04166667, 0.00000000, 0.04166667, -7.71460976, 3.86361199, 3.94921615, 3.94921884, 6.54685076],
                            [0.02747716, 0.08333333, 0.00000000, 0.08333333, -7.62519023, 3.35857616, 3.65323831, 3.65323944, 6.34969314],
                            [0.04121574, 0.12500000, 0.00000000, 0.12500000, -7.47653932, 2.71474520, 3.28609405, 3.28609571, 6.07653663]])

        self.assertTrue(np.allclose(bs._bandstructure, test_bs))

    def test_efermi_path(self):
        '''Tests the setting of the efermi filepath'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        filepath = './ut_data/total-dos-efermi.dat'

        self.assertEqual(bs._efermi_path, filepath)

    def test_kpoints(self):
        '''Tests the setting of the kpoints'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        test_kx = np.array([0.00000000, 0.04166667, 0.08333333, 0.12500000])
        test_ky = np.array([0.00000000, 0.00000000, 0.00000000, 0.00000000])
        test_kz = np.array([0.00000000, 0.04166667, 0.08333333, 0.12500000])

        kx, ky, kz = zip(*bs.kpoints)

        self.assertTrue(np.allclose(kx, test_kx))
        self.assertTrue(np.allclose(ky, test_ky))
        self.assertTrue(np.allclose(kz, test_kz))

    def test_num_bands(self):
        '''Tests the setting of the number of bands attribute'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        num_bands = 5

        self.assertEqual(bs.num_bands, num_bands)

    def test_num_kpoints(self):
        '''Tests the setting of the number of k-points attribute'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        num_kpoints = 4

        self.assertEqual(bs.num_kpoints, num_kpoints)

    def test_efermi(self):
        '''Tests the loading of the fermi energy'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        efermi = 4.070390
        self.assertEqual(bs.efermi, efermi)

    def test_energies(self):
        '''Tests the loading of the energies from the bandstructure'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)

        energies = np.array([[-7.74445402, 4.07037476, 4.07037593, 4.07038233, 6.61983340],
                             [-7.71460976, 3.86361199, 3.94921615, 3.94921884, 6.54685076],
                             [-7.62519023, 3.35857616, 3.65323831, 3.65323944, 6.34969314],
                             [-7.47653932, 2.71474520, 3.28609405, 3.28609571, 6.07653663]])

        energies = energies - 4.070390

        bs_e, bs_o = bs.get_eigenvalues()
        self.assertTrue(np.allclose(bs_e, energies))

    def test_occupancies(self):
        '''Tests the setting of the occupancies'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)
        occupancies =  np.array([[2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0],
                                [2.0, 2.0, 2.0, 2.0, 0.0]])
        bs_e, bs_o = bs.get_eigenvalues()
        self.assertTrue(np.allclose(bs_o, occupancies))

    def test_occupied_bands(self):
        '''Tests the setting of the occupied bands from the _split_bands method'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)

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

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)

        unoccupied_bands = np.array([[6.61983340, 6.54685076, 6.34969314, 6.07653663]])
        unoccupied_bands = unoccupied_bands - 4.070390
        bs_e, bs_o = bs.get_eigenvalues()
        oc, unoc = bs._split_bands()
        self.assertTrue(np.allclose(unoc, unoccupied_bands))

    def test_conduction_band(self):
        '''Tests the setting of the conduction band'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)

        # there is only one unoccupied band which also makes it the valence band
        unoccupied_bands = np.array([[6.61983340, 6.54685076, 6.34969314, 6.07653663]])
        band_min = 6.07653663

        cb, cb_min = bs._get_conduction_band(unoccupied_bands)
        self.assertTrue(np.allclose(cb, unoccupied_bands))
        self.assertEqual(cb_min, band_min)

    def test_valence_band(self):
        '''Tests the setting of the valence band'''

        bs = bandstructure.Bandstructure('Si_03082021', './ut_data/', 1.0)

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
