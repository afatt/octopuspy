#!/usr/bin/env python

'''
Generates a PROCAR and OUTCAR file using Octopus info, results, eigenvalues,
and bandstructure files

Octopus files and what needed information they contain
---------------------------------------
info -> reciprocal lattice vector, number of ions
results -> kpoint wieght
eigenvalues -> energies, occupancies, number of bands
bandstructure -> number of kpoints, num bands, CBM, VBM and  condcution, valence
bands minimum and maximum
'''

import re
import numpy as np
import numpy.ma as ma

# FILEPATH = './PbS/bands-static/'
FILEPATH = './Si/'
# FILEPATH = './New_Si/bands-static/'
# FILEPATH = './GaAs_LDA/'
INFO = FILEPATH + 'info'
RESULTS = FILEPATH + 'results.out'
BANDSTRUCTURE = FILEPATH + 'bandstructure'
EIGENVALUES = FILEPATH + 'eigenvalues'

ENERGY_SCALE = 1.0

def _calc_CBM(bandstructure):
    '''
    '''
    #_check_partial_occupancy(occupancy)
    #energy_unoccupied = ma.masked_where(occupancy > 0.5, energy)
    #return np.amin(energy_unoccupied)

    return(CBM)

def _calc_VBM(bandstructure):
    '''
    '''

    # _check_partial_occupancy(occupancy)
    # energy_occupied = ma.masked_where(occupancy < 0.5, energy)
    # return np.amax(energy_occupied)

    return(VBM)

def get_eigenvalues(bandstructure):
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
    energies = bandstructure[:, 4:]
    efermi = np.loadtxt(FILEPATH + 'total-dos-efermi.dat')[0,0]
    energies = energies - efermi

    # create a numpy array with same shape as energies
    # fill the first bands up to the valence band with an occupancy of 2.0
    occupancies = np.zeros(energies.shape)
    occupancies[energies < 0.0] = 2.0

    return(energies, occupancies)

def get_lattice_vectors(info):
    '''
    Loads the lattice vector and reciprocal lattice vector from the info file

    Params:
      info (list string): list of every line from the info
    Returns:
      zipped_vectors (zipped vectors): lattice_vector: list of length 3
                                       reciprocal_lattice_vector: list of length 3
    '''

    match = '  Lattice Vectors [A]'

    # get the first index of a list with the matched text
    #index = [idx for idx, s in enumerate(info) if match in s][0]
    index = info.index(match)

    # just get the lattice vector data and not the header text
    lattice_vector_info = info[index:index + 9]

    # extract the vector and reciprocal vector and zip together
    lattice_vector = lattice_vector_info[1:4]
    reciprocal_lattice_vector = lattice_vector_info[-3:]

    # zip for easy iteration
    zipped_vectors = zip(lattice_vector, reciprocal_lattice_vector)

    return(zipped_vectors)

def get_num_bands(bandstructure):
    '''
    Gets the number of bands by looking at the shape of the bandstructure array

    Params:
      bandstructure (numpy array): with shape (kpoints, band_info)
    Returns:
      num_bands (int): number of energy bands in the bandstructure.dat
    '''

    # numpy array in shape (kpoints, band_info) with 4 extra columns
    # # coord. kx ky kz
    num_bands = bandstructure.shape[1] - 4

    return(num_bands)

def index_of_substring(input_list, substring):
    '''
    Returns the first occurance of a substring in a list

    Params:
      input_list (list): list of string items
      substring (str): string to find the index of
    Returns:
      index (int): the index of the item containing the first occurance of
                   substring (str)
      -1 (int): there was no index found
    '''

    for index, item in enumerate(input_list):
        if substring in item:
            return(index)
    return(-1)

def get_num_kpoints(bandstructure):
    '''
    Extracts the number of kpoints from the bandstructure numpy array

    Params:
      bandstructure (numpy array): with shape (kpoints, band_info)
    Returns:
      num_kpoints (int): number of kpoints
    '''

    # numpy array in shape (kpoints, band_info)
    num_kpoints = bandstructure.shape[0]

    return(num_kpoints)

def get_num_ions(info):
    '''
    Gets the number of ions from the info list

    Params:
      info (list string): list containing every line from info file
    Returns:
      num_ions (int): number of ions
    '''

    start = ' Ion                        x              y              z'

    # TODO: this might not be the first time this occurs
    end = ' ----------------------------------------------------------'
    start_idx = info.index(start)
    end_idx = info.index(end)

    # find the length of the list between the start and end text which is
    # equivalent to number of ions. convert to string so it can be saved to file
    num_ions = str(len(info[start_idx + 1:end_idx]))

    return(int(num_ions))

def get_kpoints(bandstructure):
    '''
    Extracts the kpoints from the bandstructure numpy array

    Params:
      bandstructure (numpy array): with shape (kpoints, band_info)
    Returns:
      kpoints (zipped arrays): (kx float array, ky float array, kz float array)
    '''

    kx = bandstructure[:,1]
    ky = bandstructure[:,2]
    kz = bandstructure[:,3]

    # zip for easy iteration
    kpoints = zip(kx, ky, kz)

    return(kpoints)

def get_num_spinchannels(info):
    '''
    Calculated based on ion projection table
    '''
    return(num_spinchannels)

def gen_outcar(zipped_vectors):
    '''
    Generates the VASP OUTCAR file containing the direct lattice vector and
    reciprocal lattic vector

    Params:
      zipped_vectors (zipped vectors): lattice_vector: list of length 3
                                       reciprocal_lattice_vector: list of length 3
    '''

    f = open('OUTCAR', 'w')
    direct_header = '     direct lattice vectors'
    f.write(direct_header.ljust(46) + 'reciprocal lattice vectors\n')
    for vec, recp_vect in zipped_vectors:
        f.write(vec + ' ' + recp_vect + '\n')
    f.close()

def gen_procar(num_kpoints, num_bands, num_ions, kpoints, weights, energies, occupancies):
    '''
    TODO: Test: All kpoints must match this regular expression
    k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
    '''

    f = open('PROCAR', 'w')
    f.write('PROCAR new format' + '\n')
    f.write('# of k-points: {}          '.format(num_kpoints))
    f.write('# of bands:  {}         '.format(num_bands))
    f.write('# of ions:   {}\n\n'.format(num_ions))

    kx, ky, kz = zip(*kpoints)
    kpoints_weights = zip(kx, ky, kz, weights)

    for idx, (kx, ky, kz, weight) in enumerate(kpoints_weights):

        f.write(' k-point' + str(idx + 1).rjust(5))
        f.write(' :    ')
        f.write('{:11.8f}'.format(kx))
        f.write('{:11.8f}'.format(ky))
        f.write('{:11.8f}'.format(kz))
        f.write('     weight = {:.8f}\n\n'.format(weight))

        # num of kpoints should equal number of energy items and occupancy
        # items so idx can be used
        for i, (energy, occupancy) in enumerate(zip(energies[idx,:], occupancies[idx,:])):
            f.write('band' + '{}'.rjust(4).format(i + 1))
            f.write(' # energy' + '{:14.8f}'.format(energy))
            f.write(' # occ.' + '{:12.8f}\n\n'.format(occupancy))
            f.write('ion      s      p      d    tot\n')
            for ion in range(0, num_ions):
                f.write('  {}'.format(ion + 1))
                f.write('  {:.3f}'.format(0.000))
                f.write('  {:.3f}'.format(0.000))
                f.write('  {:.3f}'.format(0.000))
                f.write('  {:.3f}\n'.format(0.000))
            f.write('tot')
            f.write('  {:.3f}'.format(0.000))
            f.write('  {:.3f}'.format(0.000))
            f.write('  {:.3f}'.format(0.000))
            f.write('  {:.3f}\n\n'.format(0.000))
    f.close()

def get_weights(results, num_kpoints):
    '''
    Gets the kpoint weights from the results list

    Params:
      results (list string):  list of lines from results.out
      num_kpoints (int): number of kpoints
    Returns:
      weights (numpy array): shape (weights, ) dtype: float64
    '''

    # find the index of the match and get numkpoints of lines results
    match = ' index |    weight    |             coordinates              |'
    start = results.index(match)
    weights_table = results[start: start + num_kpoints + 1][1:]

    # get the second column of the table and load into numpy array
    weights_list = [row.split() for row in weights_table]
    weights_full = np.array(weights_list)
    weights = weights_full[:,2]
    weights = weights.astype('float64')

    return(weights)

def load_bandstructure():
    '''
    Loads the bandstructure file into a numpy array (without header line)

    Returns:
      bandstructure (numpy array): with shape (kpoints, band_info)
    '''

    bandstructure = np.array(np.loadtxt(BANDSTRUCTURE))

    return(bandstructure)

def load_octo_info():
    '''
    Loads the data from info file and puts into list of strings
    
    Returns:
      info (list string): list containing every line from info file
    '''
    f = open(INFO, 'r')

    # read as single string separated by newline characters
    info = f.read().splitlines()
    f.close()

    return(info)

def load_results():
    '''
    Loads the results.out file as a list of lines

    Returns:
      results (list string): list of lines from results.out
    '''

    f = open(RESULTS, 'r')
    text = f.read()

    # creates a list of lines rather than a long string with newline characters
    results = text.splitlines()
    f.close()

    return(results)

def main():
    bandstructure = load_bandstructure()
    num_bands = get_num_bands(bandstructure)
    num_kpoints = get_num_kpoints(bandstructure)
    kpoints = get_kpoints(bandstructure)

    results = load_results()
    weights = get_weights(results, num_kpoints)
    energies, occupancies = get_eigenvalues(bandstructure)

    info = load_octo_info()
    zipped_vectors = get_lattice_vectors(info)
    num_ions = get_num_ions(info)
    gen_outcar(zipped_vectors)
    gen_procar(num_kpoints, num_bands, num_ions, kpoints, weights, energies, occupancies)


if __name__ == '__main__':
    main()
