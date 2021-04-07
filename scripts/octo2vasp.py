#!/usr/bin/env python

'''
Generates a PROCAR and OUTCAR file using Octopus info, results,
and bandstructure files. The PROCAR and OUTCAR files generated are not true
VASP files, they contain only the minimum required to use in the effmass package

Octopus files and what needed information they contain
---------------------------------------
info -> reciprocal lattice vector, number of ions
results -> kpoint weight
bandstructure -> number of kpoints, num bands, energies, occupancies (calculated)
conduction band, valence band, conduction band min, and valence band max

Example use:           (required arg)       (optional arg)
python octo2vasp.py --name Si_03082021 --valence_band_index 4

Outputs:
../gen_vasp/Si_03082021/PROCAR, ../gen_vasp/Si_03082021/OUTCAR, ../gen_vasp/Si_03082021/bandstructure_plot.png
'''

import os
import re
import sys
import unicodedata
import numpy as np
import numpy.ma as ma
from glob import glob
from math import pi

# Add octopuspy to the python path
path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))

# local modules
from octopuspy.info import Info
from octopuspy.results import Results
from octopuspy.bandstructure import Bandstructure



class Octo2Vasp():

    def __init__(self, name, valence_band_index):
        self.name = name
        print("Semiconductor name: " + name)
        print("Valence Band Index: " + str(valence_band_index))
        self.filepath = self.user_prompt()
        self.bs = Bandstructure(self.filepath, valence_band_index, self.name)
        self.info = Info(self.filepath)
        self.results = Results(self.filepath, self.bs.num_kpoints)


    def gen_outcar(self):
        '''
        Generates the VASP OUTCAR file containing the direct lattice vector and
        reciprocal lattic vector. Divides the reciprocal lattice vector by 2*pi
        to match VASP OUTCAR.
        '''

        zipped_vectors = self.info.get_lattice_vectors()

        par_path = os.path.abspath(os.path.join(path, os.pardir))
        fullpaths = [file for file in glob(par_path + '/**/gen_vasp', recursive=True)]

        f = open(fullpaths[0] + '/' + self.name +  '/OUTCAR', 'w')
        direct_header = 'direct lattice vectors'
        f.write('      ')
        f.write(direct_header.ljust(39) + 'reciprocal lattice vectors\n')
        for vect, recp_vect in zipped_vectors:
            f.write('   ')
            vect = [float(v) for v in vect.split()]
            for v in vect:
                f.write('{:11.9f}'.format(v).rjust(13))
            f.write('   ')

            # VASP divides the reciprocal lattice by 2*pi
            recp_vect = [float(v)/(2*pi) for v in recp_vect.split()]
            for v in recp_vect:
                f.write('{:11.9f}'.format(v).rjust(13))

            f.write('\n')
        f.close()

    def gen_procar(self):
        '''
        Generates the VASP PROCAR file
        '''

        kpoints = self.bs.kpoints
        num_kpoints = self.bs.num_kpoints
        num_bands = self.bs.num_bands
        energies, occupancies = self.bs.get_eigenvalues()
        num_ions = self.info.num_ions
        weights = self.results.weights

        par_path = os.path.abspath(os.path.join(path, os.pardir))
        fullpaths = [file for file in glob(par_path + '/**/gen_vasp', recursive=True)]

        f = open(fullpaths[0] + '/' + self.name +  '/PROCAR', 'w')
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

    def user_prompt(self):
        '''
        Prompts the user to choose which bandstructure file to use from all
        available under the current working directory
        '''
        fullpaths = [file for file in glob('../**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        for idx, path in enumerate(fullpaths):
            print('[{}] {}\n'.format(idx + 1, path))

        filepath_choice = input('Choose the number of which bandstructure file you wish to use: ')
        while True:
            try:
                filepath_choice = int(filepath_choice) - 1
                if -1 < filepath_choice <= (len(filepaths) - 1):
                    filepath = filepaths[filepath_choice]
                    return(filepath)
                    #break
                else:
                    raise ValueError('')
            except (ValueError, IndexError) as err:
                filepath_choice = input('Choice must be number between 1 and {}, choose again: '.format(len(filepaths)))

def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens Also strip leading and trailing whitespace, dashes,
    and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    return re.sub(r'[-\s]+', '-', value).strip('-_')

def main():
    argument_list = sys.argv[1:]
    args = argument_list[::2]
    values = argument_list[1::2]

    # arg default values, raise exception on required args
    name = ''
    valence_band_index = None

    if not len(args):
        raise ValueError('Must include a name using option -n or --name')

    for curr_arg, curr_value in zip(args, values):
        if curr_arg in ('-n', '--name'):
            name = curr_value

            # chang to a valid filename
            name = slugify(name)

            # create a folder where the files for this semiconductor will be
            # saved
            try:
                os.mkdir('../gen_vasp/' + name )
            except FileExistsError as err:
                raise FileExistsError('Can not make a folder of the same name: %s' % name)

        elif curr_arg in ('-v', '--valence_band_index'):
            try:
                valence_band_index = int(curr_value)
            except ValueError as err:
                raise ValueError('Valence band index must be an integer')

    if name == '':
        raise ValueError('Must include a name using option -n or --name')

    octo2vasp = Octo2Vasp(name=name, valence_band_index=valence_band_index)
    octo2vasp.gen_outcar()
    octo2vasp.gen_procar()
    octo2vasp.bs.plot_bands()
    # octo2vasp.bs.simple_plot()


if __name__ == '__main__':
    main()
