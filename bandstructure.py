#!/usr/bin/env python

'''
'''

import numpy as np

class Bandstructure():
    ''' '''
    def __init__(self, filepath):
        self.bands = None
        self._bandstructure = None
        self._bandstructure_path = filepath + 'bandstructure'
        self._efermi_path = filepath + 'total-dos-efermi.dat'
        self.energies = None
        self.eigenvalues = None
        self.kpoints = None
        self.num_bands = 0
        self.num_kpoints = 0
        self.occupancies = None
        self._load_bandstructure()
        self._set_num_bands()
        self._set_num_kpoints()

    def get_eigenvalues(self):
        '''
        Isolates the bands from the bandstructure numpy array, subtracts the fermi
        energy and creates the energies and occupancies numpy arrays

        Params:
          bandstructure (numpy array): with shape (kpoints, band_info)
        Returns:
          energies (numpy array): with shape (num kpoints, num bands)
          occupancies (numpy array): with shape (num kpoints, num bands)
        '''

        # numpy array with shape (num kpoints, num bands)
        energies = self._bandstructure[:, 4:]
        efermi = np.loadtxt(self._efermi_path)[0,0]
        self.energies = energies - efermi

        # create a numpy array with same shape as energies
        # fill the first bands up to the valence band with an occupancy of 2.0
        self.occupancies = np.zeros(energies.shape)
        self.occupancies[energies < 0.0] = 2.0

        return(self.energies, self.occupancies)

    def get_num_bands(self):
        '''
        Gets the number of bands by looking at the shape of the bandstructure array

        Params:
          bandstructure (numpy array): with shape (kpoints, band_info)
        Returns:
          num_bands (int): number of energy bands in the bandstructure.dat
        '''

        num_bands = self.num_bands

        return(num_bands)

    def get_num_kpoints(self):
        '''
        Extracts the number of kpoints from the bandstructure numpy array

        Params:
          bandstructure (numpy array): with shape (kpoints, band_info)
        Returns:
          num_kpoints (int): number of kpoints
        '''

        num_kpoints = self.num_kpoints

        return(num_kpoints)

    def plot_bands():
        '''
        '''
        
    def _load_bandstructure(self):
        '''
        Loads the bandstructure file into a numpy array (without header line)

        Returns:
          bandstructure (numpy array): with shape (kpoints, band_info)
        '''

        self._bandstructure = np.array(np.loadtxt(self._bandstructure_path))

    def _set_num_bands(self):
        '''
        '''

        # numpy array in shape (kpoints, band_info) with 4 extra columns
        # # coord. kx ky kz
        self.num_bands = self._bandstructure.shape[1] - 4

    def _set_num_kpoints(self):
        '''
        '''

        # numpy array in shape (kpoints, band_info)
        self.num_kpoints = self._bandstructure.shape[0]
