#!/usr/bin/env python

''' Generates a PROCAR and OUTCAR file using Octopus info file'''

import re
import numpy as np
import numpy.ma as ma

# All that is used from vasppy outcar is
#  reciprocal_lattice = outcar.reciprocal_lattice_from_outcar(outcar_path)

# WHAT is used from PROCAR
# self.spin_channels = vasp_data.spin_channels
# self.number_of_bands = vasp_data.number_of_bands
# self.number_of_ions = vasp_data.number_of_ions
# number_of_kpoints = vasp_data.number_of_k_points
# vasp_data_energies = np.array( [ band.energy for band in np.ravel( vasp_data.bands ) ] )
# vasp_data_occupancies = np.array( [ band.occupancy for band in np.ravel( vasp_data.bands ) ] )

# the effective mass calculations only seem to use cartesian kpoints, and eigen values

# Files and what information they contain
# ---------------------------------------
# info -> reciprocal lattice vector, number of ions, reduced kpoints,
#         reduced kpoint weight
# results -> kpoint wieght
# eigenvalues -> energies, occupancies, number of bands
# bandstructure -> number of kpoints, num bands, CBM, VBM and  condcution, valence bands minimum and maximum

RESULTS = './PbS/bands-static/results.out'
BANDSTRUCTURE = './PbS/bands-static/bandstructure'

def get_eigen_table(info):
    '''
    '''

    start = ' #st  Spin   Eigenvalue      Occupation'
    end = 'Energy [eV]:'
    #end = '^\s*$'
    start_idx = info.index(start)
    end_idx = info.index(end)
    eigen_table = info[start_idx + 1:end_idx - 1]
    return(eigen_table)

def get_eigen_values(eigen_table):
    '''

    '''

    # remove the kpoint line info and create a list of lists consisting of
    # each row and column of the eigen table
    k_removed_list = [line for line in eigen_table if '#k' not in line]
    eigen_list_split = [line.split() for line in k_removed_list]

    # save the engergies and occupancies to a new list if the line is not empty
    energies = [float(line[2]) for line in eigen_list_split if line]
    occupancies = [float(line[3]) for line in eigen_list_split if line]

    return(energies, occupancies)

def group_eigen_values(num_bands, energies, occupancies):
    '''
        TODO: transfer to ndarray with shape (num kpoints, num bands)
    '''

    # array of form [[kpoint1_energies], ... [kpoint16_energies]]
    # array of form [[kpoint1_occupancies], ... [kpoint16_occupancies]]
    grouped_energies = [energies[i:i + num_bands] for i in range(0, len(energies), num_bands)]
    grouped_occupancies = [occupancies[i:i + num_bands] for i in range(0, len(occupancies), num_bands)]

    return(grouped_energies, grouped_occupancies)

def eigen_values_np(grouped_energies, grouped_occupancies):
    '''TODO: lump the three functions to work together and use np arrays in main code
    '''

    # numpy array with shape (num bands, num kpoints)
    energies = np.array(grouped_energies).T
    occupancies = np.array(grouped_occupancies).T

    # mask in energies where occupancies is > 0.5
    energy_unoccupied = ma.masked_where(occupancies > 0.5, energies)
    print(energy_unoccupied)
    print(np.amax(energy_unoccupied))

    energy_occupied = ma.masked_where(occupancies < 0.5, energies)
    print(np.amax(energy_occupied))

    return(energies, occupancies)

def get_lattice_vectors(info):
    '''TODO: Original lattice vector had 8 decimal places
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
    '''
    start = ' Ion                        x              y              z'

    # TODO: this might not be the first time this occurs
    end = ' ----------------------------------------------------------'
    start_idx = info.index(start)
    end_idx = info.index(end)

    # find the length of the list between the start and end text which is
    # equivalent to number of ions. convert to string so it can be saved to file
    num_ions = str(len(info[start_idx + 1:end_idx]))
    return(num_ions)

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
    f = open('OUTCAR', 'w')
    direct_header = '     direct lattice vectors'
    f.write(direct_header.ljust(46) + 'reciprocal lattice vectors\n')
    for vec, recp_vect in zipped_vectors:
        f.write(vec + ' ' + recp_vect + '\n')
    f.close()

def gen_procar(num_kpoints, num_bands, num_ions, kpoints, energies, occupancies):
    '''
    TODO: Test: All kpoints must match this regular expression
    k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
    '''

    f = open('PROCAR', 'w')
    f.write('PROCAR new format' + '\n')
    f.write('# of k-points: {}          # of bands:  {}         # of ions:   {}\n\n'.format(num_kpoints, num_bands, num_ions))

    e_group, o_group = group_eigen_values(num_bands, energies, occupancies)

    x, y = eigen_values_np(e_group, o_group)

    for idx, (kx, ky, kz, weight) in enumerate(kpoints): # 16 points

        f.write(' k-point' + str(idx + 1).rjust(5))
        f.write(' :    ')
        f.write('{:11.8f}'.format(float(kx)))
        f.write('{:11.8f}'.format(float(ky)))
        f.write('{:11.8f}'.format(float(kz)))
        f.write('     weight = {:.8f}\n\n'.format(float(weight)))

        # num of kpoints should equal number of energy groups and occupancy
        # groups so idx can be used
        for i, (energy, occupancy) in enumerate(zip(e_group[idx], o_group[idx])):
            f.write('band' + '{}'.rjust(4).format(i + 1))
            f.write(' # energy' + '{:14.8f}'.format(energy))
            f.write(' # occ.' + '{:12.8f}\n\n'.format(occupancy))
            f.write('ion      s      p      d    tot\n')
            f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(1, 0.065, 0.000, 0.000, 0.065))
            f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(2, 0.556, 0.000, 0.000, 0.556))
            f.write('tot  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n\n'.format(0.620, 0.000, 0.000, 0.620))
    f.close()

    # wf-kpoint#-band#.cube
    # wf-k001-st0002.cube
    # orbitals = [[s, p, d], [s, py, pz, px, dxy, dyz, dz2, dxz, dx2]]

def get_weights(results, num_kpoints):
    '''
    Gets the kpoint weights from the results list

    Params:
      results (list string):  list of lines from results.out
      num_kpoints (int): number of kpoints
    Returns:
      weights (numpy array): shape (weights, )
    '''

    # find the index of the match and get numkpoints of lines results
    match = ' index |    weight    |             coordinates              |'
    start = results.index(match)
    weights_table = results[start: start + num_kpoints][1:]

    # get the second column of the table and load into numpy array
    weights_list = [row.split() for row in weights_table]
    weights_full = np.array(weights_list)
    weights = weights_full[:,2]

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
    '''
    f = open('./PbS/bands-static/info', 'r')

    # read as single string separated by newline characters
    info = f.read()

    # creates a list of lines rather than a long string with newline characters
    info_list = info.splitlines()
    f.close()

    return(info_list, info)

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

    # info_list, info = load_octo_info()

    # zipped_vectors = get_lattice_vectors(info_list)
    #
    # num_ions = get_num_ions(info_list)
    # eigen_table = get_eigen_table(info_list)
    # energies, occupancies = get_eigen_values(eigen_table)
    # gen_outcar(zipped_vectors)
    # gen_procar(num_kpoints, num_bands, num_ions, kpoints, energies, occupancies)


if __name__ == '__main__':
    main()
