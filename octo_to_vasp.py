#!/usr/bin/env python

''' Generates a PROCAR and OUTCAR file using Octopus info file'''

import re
import numpy as np

# All that is used from vasppy outcar is
#  reciprocal_lattice = outcar.reciprocal_lattice_from_outcar(outcar_path)

# WHAT is used from PROCAR
# self.spin_channels = vasp_data.spin_channels
# self.number_of_bands = vasp_data.number_of_bands
# self.number_of_ions = vasp_data.number_of_ions
# number_of_kpoints = vasp_data.number_of_k_points
# vasp_data_energies = np.array( [ band.energy for band in np.ravel( vasp_data.bands ) ] )
# vasp_data_occupancies = np.array( [ band.occupancy for band in np.ravel( vasp_data.bands ) ] )

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
    energies = [line[2] for line in eigen_list_split if line]
    occupancies = [line[3] for line in eigen_list_split if line]

    return(energies, occupancies)

def group_eigen_values(num_bands, energies, occupancies):
    '''
    '''

    grouped_energies = [energies[i:i + num_bands] for i in range(0, len(energies), num_bands)]
    grouped_occupancies = [occupancies[i:i + num_bands] for i in range(0, len(occupancies), num_bands)]

    # [([energy_list1], [occ_list1]), ([energy_list2], [occ_list2])]
    #eigen_groups = zip(grouped_energies, grouped_occupancies)

    return(grouped_energies, grouped_occupancies)

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

def get_num_bands(info):
    '''
    Gets the number of bands by looking at the eigen value table in info

    Returns:
      num_bands (int): number of energy bands in the bandstructure.dat
    '''

    eigen_table = get_eigen_table(info)
    start_idx = index_of_substring(eigen_table, '#k =   1')
    end_idx = index_of_substring(eigen_table, '#k =   2')
    num_bands = (end_idx - start_idx) - 1

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

def get_num_kpoints(info):
    '''
    '''

    # TODO: Move to a separate function
    match = 'Number of symmetry-reduced k-points'
    matched = [x for x in info if match in x]

    # parse into list with space delimiter and get the last item
    kpoints_line = matched[0].split()
    num_kpoints = kpoints_line[-1]

    #TODO: Check to make sure its an integer value

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

def get_kpoints(info, num_kpoints):
    '''TODO: The data from eigen values has more decimal points
    '''

    # get the kpoints data from the info data
    # there are three text header lines that need to be skipped
    kpoint_idx = info.index('List of k-points:') + 3
    kpoint_info = info[kpoint_idx:kpoint_idx + int(num_kpoints)]

    # split each item in list with space delimiter and remove the kpoint number
    # of each item
    kpoints = [x.split()[1:] for x in kpoint_info]

    # create individual lists of each kpoint item
    kx  = [kp[0] for kp in kpoints]
    ky = [kp[1] for kp in kpoints]
    kz = [kp[2] for kp in kpoints]
    weight = [kp[3] for kp in kpoints]

    # zip for easy iteration
    kpoints = zip(kx, ky, kz, weight)

    return(kpoints)

def get_num_spinchannels(info):
    '''
    Calculated based on ion projection table
    '''
    return(num_spinchannels)

def load_band_data():
    '''
    '''
    data = np.array(np.loadtxt('./PbS/bands-static/bandstructure'))
    return(data)

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
            f.write(' # energy' + '{:14.8f}'.format(float(energy)))
            f.write(' # occ.' + '{:12.8f}\n\n'.format(float(occupancy)))
            f.write('ion      s      p      d    tot\n')
            f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(1, 0.065, 0.000, 0.000, 0.065))
            f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(2, 0.556, 0.000, 0.000, 0.556))
            f.write('tot  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n\n'.format(0.620, 0.000, 0.000, 0.620))
    f.close()

    # wf-kpoint#-band#.cube
    # wf-k001-st0002.cube
    # orbitals = [[s, p, d], [s, py, pz, px, dxy, dyz, dz2, dxz, dx2]]


def main():
    info_list, info = load_octo_info()
    num_bands = get_num_bands(info_list)
    zipped_vectors = get_lattice_vectors(info_list)
    num_kpoints = get_num_kpoints(info_list)
    num_ions = get_num_ions(info_list)
    kpoints = get_kpoints(info_list, num_kpoints)
    eigen_table = get_eigen_table(info_list)
    energies, occupancies = get_eigen_values(eigen_table)
    gen_outcar(zipped_vectors)
    gen_procar(num_kpoints, num_bands, num_ions, kpoints, energies, occupancies)


if __name__ == '__main__':
    main()
