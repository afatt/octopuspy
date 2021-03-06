'''
'''

import numpy as np
import matplotlib.pyplot as plt

class Bandstructure():
    ''' '''
    def __init__(self, filepath, energy_scale):
        self._bandstructure = None
        self._bandstructure_path = filepath + 'bandstructure'
        self._efermi_path = filepath + 'total-dos-efermi.dat'
        self.efermi = np.loadtxt(self._efermi_path)[0,0]
        self.eigenvalues = None
        self.energies = None
        self.energy_scale = energy_scale
        self.kpoints = None
        self.num_bands = 0
        self.num_kpoints = 0
        self.occupancies = None
        self._load_bandstructure()
        self._set_num_bands()
        self._set_num_kpoints()
        self._set_kpoints()


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
        energies = energies - self.efermi
        self.energies = energies / self.energy_scale

        # create a numpy array with same shape as energies
        # fill the first bands up to the valence band with an occupancy of 2.0
        self.occupancies = np.zeros(energies.shape)
        self.occupancies[energies < 0.0] = 2.0

        return(self.energies, self.occupancies)

    def get_kpoints(self):
        '''
        Extracts the kpoints from the bandstructure numpy array

        Params:
          bandstructure (numpy array): with shape (kpoints, band_info)
        Returns:
          kpoints (zipped arrays): (kx float array, ky float array, kz float array)
        '''

        kpoints = self.kpoints

        return(kpoints)

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

    def plot_bands(self):
        '''
        Creates and saves a figure of the bandstructure plots. The valence band
        is colored darkorange and the conduction band is colored yellow.The
        valence band max and conduction band minimum are labeled with a '*'

        Output:
          bandstructure_plot.png: bandstructure figure saved to the current
                                  working directory
        '''

        x_data = self._bandstructure[:,0]

        occupied_bands, unoccupied_bands = self._split_bands()
        valence_band, vb_max = self._get_valence_band(occupied_bands)
        conduction_band, cb_min= self._get_conduction_band(unoccupied_bands)
        vb_max_index = valence_band.argmax()
        cb_min_index = conduction_band.argmin()

        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])

        # create plot for all but last row (valence band)
        # plot valence band on its own in color darkorange
        for band in occupied_bands[:-1,:]:
            ax.plot(x_data, band, color='darkred')
        ax.plot(x_data, valence_band, color='darkorange')

        # create plot for all but first row (conduction band)
        # plot conduction band on its own in yellow
        for band in unoccupied_bands[1:,:]:
            ax.plot(x_data, band, color='indigo')
        ax.plot(x_data, conduction_band, color='yellow')

        ax.plot(x_data[vb_max_index], vb_max, '*')
        ax.plot(x_data[cb_min_index], cb_min, '*')
        ax.set_ylabel('E-Ef')
        ax.set_title('Bulk Bandstructure')
        ax.set_xticks([0.00,0.118,0.255,0.391,0.726,1.000])
        ax.set_xticklabels(['K','Gamma','X','W','K','Gamma','L','U','W','L','K'])
        ax.set_xticks([0.00,0.148,0.289,0.359,0.408,0.557,0.679,0.765,0.814,0.914,1.000])
        ax.tick_params(axis='both',labelsize=12)

        plt.axhline(y=0)
        fig.savefig('bandstructure_plot.png')

    def _get_conduction_band(self, unoccupied_bands):
        '''
        Separates the conduction band into its own numpy array

        Args:
          unoccupied_bands (numpy array): shape (num_bands, num_kpoints)
        Returns:
          conduction_band (numpy array): shape (num_kpoints, )
          cb_min (float): The smallest value in the conduction_band numpy array
        '''

        conduction_band = unoccupied_bands[0, :]
        cb_min = conduction_band.min()
        return(conduction_band, cb_min)

    def _get_valence_band(self, occupied_bands):
        '''
        Separates the valence band into its own numpy array

        Args:
          occupied_bands (numpy array): shape (num_bands, num_kpoints)
        Returns:
          valence_band (numpy array): shape (num_kpoints, )
          vb_max (float): The largest value in the valence_band numpy array
        '''

        valence_band = occupied_bands[-1, :]
        vb_max = valence_band.max()
        return(valence_band, vb_max)

    def _load_bandstructure(self):
        '''
        Loads the bandstructure file into a numpy array (without header line)

        Returns:
          bandstructure (numpy array): with shape (kpoints, band_info)
        '''

        self._bandstructure = np.array(np.loadtxt(self._bandstructure_path))

    def _set_kpoints(self):
        '''
        '''

        kx = self._bandstructure[:,1]
        ky = self._bandstructure[:,2]
        kz = self._bandstructure[:,3]

        # zip for easy iteration
        self.kpoints = zip(kx, ky, kz)

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

    def _split_bands(self):
        '''
        Takes the energies numpy array shape (num_bands, num_kpoints) and
        creates two numpy arrays, occupied and unoccupied bands

        Returns:
          occupied_bands (numpy array): shape (num_bands, num_kpoints)
          unoccupied_bands (numpy array): shape (num_bands, num_kpoints)
        '''

        # transpose to shape (num_kpoints, num_bands) before separating into
        # occupied and unoccupied bands of shape (num_energies, )
        occupied_bands = self.energies.T[self.energies.T < 0.0]
        unoccupied_bands = self.energies.T[self.energies.T >= 0.0]

        # the number of occupied/unoccupied bands can be found from the shape
        # of those numpy arrays
        num_occupied_bands = int(occupied_bands.shape[0] / self.num_kpoints)
        num_unoccupied_bands = int(unoccupied_bands.shape[0] / self.num_kpoints)

        # reshape so its not just a single dimensional numpy array, and it has
        # a row for each band
        occupied_bands = np.reshape(occupied_bands,
                                   (num_occupied_bands, self.num_kpoints))

        unoccupied_bands = np.reshape(unoccupied_bands,
                                     (num_unoccupied_bands, self.num_kpoints))

        return(occupied_bands, unoccupied_bands)
