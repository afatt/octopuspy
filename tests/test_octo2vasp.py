import os
import re
import sys
import unittest
import numpy as np
from glob import glob
from io import StringIO

from scripts.octo2vasp import Octo2Vasp

test_data_dir = 'test_data'
if 'win' in sys.platform:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '\\'
else:
    test_path = os.path.join( os.path.dirname( __file__ ), test_data_dir ) + '/'

class TestOcto2Vasp(unittest.TestCase):
    ''' '''

    def test_gen_procar_pos_kpoints(self):
        '''
        Test Case 1: positive kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('../**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '..\\octopuspy\\tests\\test_data\\pos_kpoints\\bandstructure'
        lin_path = '../octopuspy/tests/test_data/pos_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', occ_band_num=4)

        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')

        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                # print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)

    def test_gen_procar_neg_kpoints(self):
        '''
        Test Case 2: negative kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('../**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '..\\octopuspy\\tests\\test_data\\neg_kpoints\\bandstructure'
        lin_path = '../octopuspy/tests/test_data/neg_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', occ_band_num=4)

        o2v.filepath = './test_data'
        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')
        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                # print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)

    def test_gen_procar_mixed_sign_kpoints(self):
        '''
        Test Case 3: mixed sign kpoints
        All kpoints must match this regular expression
        k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)
        '''

        fullpaths = [file for file in glob('../**/bandstructure*', recursive=True)]
        filepaths = [os.path.dirname(path) + '/' for path in fullpaths]
        win_path = '..\\octopuspy\\tests\\test_data\\mixed_sign_kpoints\\bandstructure'
        lin_path = '../octopuspy/tests/test_data/mixed_sign_kpoints/bandstructure'
        try:
            idx = fullpaths.index(win_path)
        except Exception as err:
            print('Could not match windows path trying Linux')

        try:
            idx = fullpaths.index(lin_path)
        except Exception as err:
            print('Could not match linux path, windows will be used')

        selection = idx + 1

        oldstdin = sys.stdin
        sys.stdin = StringIO(str(selection))
        o2v = Octo2Vasp(name='test_out', occ_band_num=4)

        o2v.filepath = './test_data'
        try:
            os.mkdir('./gen_vasp/test_out')
        except FileExistsError as err:
            print('./gen_vasp/test_out already exists')
        o2v.gen_procar()

        # the number of lines that should match the kpoint regular expression
        num_matched = 4
        pattern = re.compile('k-point\s+(\d+)\s*:\s+([- ][01].\d{8})([- ][01].\d{8})([- ][01].\d{8})\s+weight = ([01].\d+)')

        matched = []
        f = open('./gen_vasp/test_out/PROCAR', 'r')
        for i, line in enumerate(f):
            for match in re.finditer(pattern, line):
                print('Found on line %s: %s' % (i+1, match.group()))
                matched.append(match.group)

        # close the filehandle
        f.close()
        os.remove('./gen_vasp/test_out/PROCAR')
        os.rmdir('./gen_vasp/test_out')

        self.assertEqual(len(matched), num_matched)


if __name__ == '__main__':
    unittest.main()
