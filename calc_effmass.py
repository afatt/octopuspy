#!/usr/bin/env python

'''
Using the example code from the Jupyter notebook:
https://nbviewer.jupyter.org/github/lucydot/effmass/blob/master/paper/notebook.ipynb
This script is meant to provide the results of the effective mass of any semiconductor
'''

import sys
from shutil import copyfile

# import scientific libraries
import math
import matplotlib.pyplot as plt
import numpy as np

from effmass import extrema

# import modules from the effmass package
from effmass import inputs, analysis, extrema, outputs, dos, ev_to_hartree

extrema_search_depth = 0.075
energy_range = 0.75
ignored_kpoints = 0

semicond = input('type the name of your semiconductor: ')

#settings = inputs.Settings(extrema_search_depth=0.075, energy_range=0.25)
settings = inputs.Settings(extrema_search_depth=extrema_search_depth,
                           energy_range=energy_range)

# PROCAR = "./data/CdTe/HSE06_SoC/PROCAR"
# PROCAR = "./data/CdTe/HSE06/PROCAR"
# PROCAR = "./data/GaAs/HSE06/PROCAR"
# PROCAR = "./data/GaAs/LDA/PROCAR"
# PROCAR = "./data/MAPI/HSE06/PROCAR"
PROCAR = "PROCAR"

#data = inputs.Data("./data/CdTe/HSE06_SoC/OUTCAR", PROCAR, ignore=0)
#data = inputs.Data("./data/CdTe/HSE06/OUTCAR",, ignore=0)
#data = inputs.Data("./data/GaAs/HSE06/OUTCAR", PROCAR, ignore=0)
#data = inputs.Data("./data/GaAs/LDA/OUTCAR",PROCAR, ignore=0)
#data = inputs.Data("./data/MAPI/HSE06/OUTCAR", PROCAR, ignore=0)
data = inputs.Data("OUTCAR", PROCAR, ignore=ignored_kpoints)

ans = extrema.find_extrema_indices(data, settings)
# print(ans)
# for band, kpoint in ans:
#     print('band: ' + str(band))
#     print('kpoint: ' + str(kpoint))


print('VBM: ' + str(data.VBM))
print('CBM: ' + str(data.CBM))
print('Band Gap: ' + str(data.CBM - data.VBM))

print('Fermi Energy: ' + str(data.fermi_energy))
print('num bands: ' + str(data.number_of_bands))
print('num ions: ' + str(data.number_of_ions))
print('len_kpoints: ' + str(len(data.kpoints)))
print('len_energies: ' + str(len(data.energies)))
print('len_occ: ' + str(len(data.occupancy)))
print(data.energies)
print(data.occupancy)

segments = extrema.generate_segments(settings,data)
print(segments)
#str(segments[-1])

fig, ax = outputs.plot_segments(data,settings,segments)
fig.savefig('./results/' + semicond + '_segments1.png')

#fig, ax = outputs.plot_segments(data,settings,[segments[-1],segments[-3]])
#fig.savefig('./results/' + semicond + '_segments2.png')

for segment in segments:
    try:
        print('Effmass: ' + str(segment.five_point_leastsq_effmass()))
    except Exception as err:
        print(err)

f = open('./results/effmass_of_' + semicond + '.txt', 'w')
leastsq_effmass = segments[-1].five_point_leastsq_effmass()
f.write('3-point parabolic mass: {}\n'.format(leastsq_effmass))
print(leastsq_effmass)

finite_diff_effmass = segments[-1].finite_difference_effmass()
f.write('3-point finite difference mass: {}\n'.format(finite_diff_effmass))
print(finite_diff_effmass)

parabolic_effmass = segments[-1].weighted_leastsq_effmass()
f.write('weighted parabolic mass: {}\n'.format(parabolic_effmass))
print(parabolic_effmass)
f.close()

f = open('./results/summary_' + semicond + '.txt', 'w')
f.write('ignored kpoints used: ' + str(ignored_kpoints) + '\n')
f.write('extrema search depth used: ' + str(extrema_search_depth) + '\n')
f.write('energy range used: ' + str(energy_range) + '\n')
f.write('VBM: ' + str(data.VBM) + '\n')
f.write('CBM: ' + str(data.CBM) + '\n')
f.write('Bandgap: ' + str(data.CBM - data.VBM) + '\n')
f.close()

original_stdout = sys.stdout # Save a reference to the original standard output

with open('./results/summary_' + semicond + '.txt', 'a') as f:
    sys.stdout = f # Change the standard output to the file we created.
    outputs.print_results(segments[-1], data, settings)
    sys.stdout = original_stdout # Reset the standard output to its original value

# make a copy of the PROCAR File used
copyfile(PROCAR, './results/PROCAR')
