# octopuspy

[![Build Status](https://travis-ci.com/afatt/octopuspy.svg?token=vuSwZWN5GqEyqapXM9jc&branch=master)](https://travis-ci.com/afatt/octopuspy)

[![Test Coverage](https://api.codeclimate.com/v1/badges/4fc9872ae9f9c1327c21/test_coverage)](https://codeclimate.com/github/afatt/octopuspy/test_coverage)

`octopuspy` is a Python package that allows for the manipulation and processing of files generated by the [Octopus](https://octopus-code.org/wiki/Main_Page) DFT scientific program. This package was created with the
intent of using with the [effmass](https://github.com/lucydot/effmass) package.

## Installation

`octopuspy` can be installed from PyPI:

    pip install octopuspy

or if the package exists locally:

    cd octopuspy
    python setup.py install

## How To Use

Required Files: bandstructure, info, results.out, and dos-xxxx.dat files (if these don't exists the number of occupied bands can be specified)

**bandstructure**
file includes: number of kpoints, number of bands, energies, and occupancies (calculated). The Bandstructure class uses the dos-xxx.dat files to determine the number of occupied bands, if none are found user will be prompted to use the `occ_band_num` to set the number of occupied bands:

    band_data = Bandstructure(name, bandstructure_filepath, occ_band_num=10)

**info** file includes: direct lattice vector, reciprocal lattice vector, and number of ions:

    info_data = Info(info_filepath)

**results.out** file includes: kpoint weights

    results_out = Results(results_filepath, num_kpoints)

The `octopuspy` package includes the script `octo2vasp.py` for converting the  

### Recommended Octopus Settings

## Testing

Unittests are included with this package and can be run using:

    python -m unittest discover
