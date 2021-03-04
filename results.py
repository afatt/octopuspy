'''
'''

import numpy as np

class Results():
    ''' '''
    def __init__(self, filepath, num_kpoints):
        self.weights = None
        self._num_kpoints = num_kpoints
        self._results = None
        self._results_path = filepath + 'results.out'
        self._load_results()
        self.get_weights()

    def get_weights(self):
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
        start = self._results.index(match)
        weights_table = self._results[start: start + self._num_kpoints + 1][1:]

        # get the second column of the table and load into numpy array
        weights_list = [row.split() for row in weights_table]
        weights_full = np.array(weights_list)
        weights = weights_full[:,2]
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
