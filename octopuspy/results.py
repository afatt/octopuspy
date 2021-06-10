'''
Uses the results.out file generated by Octopus and makes its information
available to use. results.out is not provided by Octopus by default. User
must redirect the stdout of Octopus to a new file named results.out while
executing the Octopus program. Octopus > results.out

Information contained in results.out file
-------------------------------------------
results -> kpoint weight
'''

import numpy as np

class Results():
    '''
    Class that holds and gives methods to the information of a results.out
    file

    Attribues:
      weights (numpy array): with shape (num_kpoints, ) dtype='float64'
      _num_kpoints (int): number of kpoints
      _results (list string): list of lines from results.out
      _results_path (string): filepath with the addition of the results.out file
    '''
    def __init__(self, filepath, num_kpoints):
        '''
        Args:
          filepath (string): the filepath to the results.out file
        '''
        self.weights = None
        self._num_kpoints = num_kpoints
        self._results = None
        self._results_path = filepath + 'results.out'
        self._load_results()
        self.get_weights()

    def get_weights(self):
        '''
        Gets the kpoint weights from the results list

        Returns:
          weights (numpy array): shape (weights, ) dtype: float64
        '''

        # find the index of the match and get numkpoints of lines results
        match = ' index |    weight    |             coordinates              |'
        start = self._results.index(match)

        # start at the second line not including the matched text
        weights_table = self._results[start: start + 2][1:]

        # get the second column of the table and load into numpy array
        line = [row.split() for row in weights_table]
        weight = line[0][2]

        # create a 1 dimensional numpy array with type of float
        weights = np.full((self._num_kpoints, ), weight)
        self.weights = weights.astype('float64')

        return(self.weights)

    def _load_results(self):
        '''
        Loads the results.out file as a list of lines

        Returns:
          results (list string): list of lines from results.out
        '''

        f = open(self._results_path, 'r')
        text = f.read()

        # creates a list of lines rather than a long string with newline characters
        self._results = text.splitlines()
        f.close()
