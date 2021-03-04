'''
'''

class Info():
    ''' '''
    def __init__(self, filepath):
        self._info = None
        self._info_path = filepath + 'info'
        self.lattice_vector = None
        self.reciprocal_lattice_vector = None
        self.num_ions = 0
        self._load_info()
        self.get_num_ions()


    def get_num_ions(self):
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
        start_idx = self._info.index(start)
        end_idx = self._info.index(end)

        # find the length of the list between the start and end text which is
        # equivalent to number of ions. convert to string so it can be saved to file
        num_ions = str(len(self._info[start_idx + 1:end_idx]))

        self.num_ions = int(num_ions)

        return(self.num_ions)

    def get_lattice_vectors(self):
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
        index = self._info.index(match)

        # just get the lattice vector data and not the header text
        lattice_vector_info = self._info[index:index + 9]

        # extract the vector and reciprocal vector and zip together
        self.lattice_vector = lattice_vector_info[1:4]
        self.reciprocal_lattice_vector = lattice_vector_info[-3:]

        # zip for easy iteration
        zipped_vectors = zip(self.lattice_vector, self.reciprocal_lattice_vector)

        return(zipped_vectors)

    def _load_info(self):
        '''
        Loads the data from info file and puts into list of strings

        Returns:
          info (list string): list containing every line from info file
        '''
        f = open(self._info_path, 'r')

        # read as single string separated by newline characters
        self._info = f.read().splitlines()
        f.close()
