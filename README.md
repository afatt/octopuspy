# octopuspy

[![Build Status](https://travis-ci.com/afatt/octopuspy.svg?token=vuSwZWN5GqEyqapXM9jc&branch=master)](https://travis-ci.com/afatt/octopuspy)
[![Test Coverage](https://api.codeclimate.com/v1/badges/4fc9872ae9f9c1327c21/test_coverage)](https://codeclimate.com/github/afatt/octopuspy/test_coverage)

`octopuspy` is a Python package that allows for the manipulation and processing of files generated by the [Octopus](https://octopus-code.org/wiki/Main_Page) DFT scientific program. This package was created with the intent of using with the [effmass](https://github.com/lucydot/effmass) package.

## Installation

`octopuspy` can be installed from PyPI:

    pip install octopuspy

Or download the source code and install:

    git clone https://github.com/afatt/octopuspy.git
    cd octopuspy
    python3 setup.py install

## How To Use

Required `Octopus` files: **bandstructure**, **info**, **results.out**(typically out.log), **eigenvalues**, and **total-dos-efermi.dat**

### API

**bandstructure**
file includes: number of kpoints, number of bands, energies:

    band_data = Bandstructure(filepath='./GaAs_HSE06', name='GaAs_HSE06_04092021')

**info** file includes: direct lattice vector, reciprocal lattice vector, and number of ions:

    info_data = Info(filepath='./GaAs_HSE06')

**results.out** file includes: kpoint weights. **results.out** is obtained by redirecting the stdout to a new file named results.out when running **Octopus**.`Octopus > results.out`. More explanation found on the Octopus wiki [here](https://octopus-code.org/wiki/Manual:Running_Octopus)

    results_out = Results(filepath='./GaAs_HSE06', band_data.num_kpoints)

### Scripts

The `octopuspy` package includes the script `octo2vasp.py` for converting the **bandstructure**, **info**, **results.out**, **eigenvalues** files into [VASP](https://www.vasp.at/) PROCAR and OUTCAR files.⚠️ *The PROCAR and OUTCAR files generated are not true VASP files, they contain only the minimum data required to use in the effmass package*.

Copy your `Octopus` files: **bandstructure**, **info**, **results.out**, and **eigenvalues** into a new folder anywhere inside the octopuspy package and run:

	cd scripts
	python3 octo2vasp.py --name name_of_semiconductor

The `octo2vasp.py` script will write the PROCAR, OUTCAR, and a bandstructure_plot.png to the `gen_vasp` folder inside the octopuspy package.

### Recommended Octopus Settings

Garbage in equals garbage out and testing this package produced a lot of garbage. Included with this package are Octopus settings that helped to achieve results matching the literature **anotated_example_inp**. The **inp** file is also included with each set of `/test_data` for reference. More information about input settings can be found on the Octopus wiki [here](https://octopus-code.org/wiki/Manual:Input_file)

## Testing

Unittests are included with this package and can be run using:

    python -m unittest discover
