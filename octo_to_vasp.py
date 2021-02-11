#!/usr/bin/env python

''' Generates a PROCAR and OUTCAR file using Octopus output data'''

import re
import numpy as np

# All that is used from vasppy outcar is
#  reciprocal_lattice = outcar.reciprocal_lattice_from_outcar(outcar_path)

# WHAT is used from PROCAR
# self.spin_channels = vasp_data.spin_channels
# self.number_of_bands = vasp_data.number_of_bands
# self.number_of_ions = vasp_data.number_of_ions
# number_of_kpoints = vasp_data.number_of_k_points

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

    zipped_vectors = zip(lattice_vector, reciprocal_lattice_vector)

    return(zipped_vectors)

def get_num_bands():
    '''
    Opens bandstructure.dat and reads the first line of the file and gets the
    number of bands

    Returns:
      num_bands (str): number of energy bands in the bandstructure.dat
    '''

    # read the first line that containes the nmber of bands
    f = open('bandstructure.dat', 'r')
    info_line = f.readline()
    f.close()

    # create a list out of items in the first line
    # find the index of bands: and the next index is the number of bands
    info_list = info_line.split()
    num_bands = info_list[info_list.index('bands:') + 1]
    return(num_bands)

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
    '''
    '''

    # get the kpoints data from the info data
    # there are three text header lines that need to be skipped
    kpoint_idx = info.index('List of k-points:') + 3
    kpoint_info = info[kpoint_idx:kpoint_idx + int(num_kpoints)]

    # split each item in list with space delimiter and remove the kpoint number
    # of each item
    kpoints = [x.split()[1:] for x in kpoint_info]

    # zip for easy iteration
    kpoints = zip(kpoints)

    return(kpoints)

def get_num_spinchannels(info):
    '''
    '''
    return(num_spinchannels)

def load_band_data():
    '''
    '''
    data = np.array(np.loadtxt('bandstructure.dat'))
    return(data)

def load_octo_info():
    '''
    '''
    f = open('./PbS/bands-static/info', 'r')

    # creates a list of lines rather than a long string with newline characters
    info = f.read().splitlines()
    f.close()
    return(info)

def main():
    num_bands = get_num_bands()
    info = load_octo_info()
    zipped_vectors = get_lattice_vectors(info)
    num_kpoints = get_num_kpoints(info)
    num_ions = get_num_ions(info)
    kpoints = get_kpoints(info, num_kpoints)
    print(kpoints)


    f = open('OUTCAR', 'w')
    f.write("      direct lattice vectors                 reciprocal lattice vectors\n")
    for vec, recp_vect in zipped_vectors:
        f.write(' ' + vec + ' ' + recp_vect + '\n')
    f.close()

    f = open('PROCAR', 'w')
    f.write('PROCAR new format' + '\n')
    f.write('# of k-points: {}          # of bands:  {}         # of ions:   {}\n\n'.format(num_kpoints, num_bands, num_ions))
    f.write(' k-point    1 :    {:.8f} {:.8f} {:.8f}     weight = {:.8f}\n\n'.format(0.00000000, 0.00000000, 0.00000000, 0.00462963))
    f.write(' band   {} # energy  {:.8f} # occ.  {:.8f}\n\n'.format(1, -10.30320759, 1.00000000))
    f.write('ion      s      p      d    tot\n')
    f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(1, 0.065, 0.000, 0.000, 0.065))
    f.write('  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(2, 0.556, 0.000, 0.000, 0.556))
    f.write('tot  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n'.format(0.620, 0.000, 0.000, 0.620))
    f.close()


    # k - point 1 - 465
      # loop through 1 - 48 bands
        # ions 1, 2 (four different times?)
'''
ion      s      p      d    tot
  1  0.065  0.000  0.000  0.065
  2  0.556  0.000  0.000  0.556
tot  0.620  0.000  0.000  0.620
  1  0.037 -0.000 -0.000  0.037
  2  0.320 -0.000 -0.000  0.320
tot  0.358 -0.000 -0.000  0.358
  1  0.037 -0.000  0.000  0.037
  2  0.320  0.000  0.000  0.320
tot  0.358  0.000  0.000  0.358
  1  0.037  0.000  0.000  0.037
  2  0.322  0.000  0.000  0.322
tot  0.359  0.000  0.000  0.359
'''
if __name__ == '__main__':
    main()
