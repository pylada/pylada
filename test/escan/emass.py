###############################
#  This file is part of PyLaDa.
#
#  Copyright (C) 2013 National Renewable Energy Lab
# 
#  PyLaDa is a high throughput computational platform for Physics. It aims to make it easier to submit
#  large numbers of jobs on supercomputers. It provides a python interface to physical input, such as
#  crystal structures, as well as to a number of DFT (VASP, CRYSTAL) and atomic potential programs. It
#  is able to organise and launch computational jobs on PBS and SLURM.
# 
#  PyLaDa is free software: you can redistribute it and/or modify it under the terms of the GNU General
#  Public License as published by the Free Software Foundation, either version 3 of the License, or (at
#  your option) any later version.
# 
#  PyLaDa is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
#  Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along with PyLaDa.  If not, see
#  <http://www.gnu.org/licenses/>.
###############################

from math import fabs as abs
from numpy import array
from pylada.crystal import fill_structure
from pylada.escan import read_input
from pylada.escan.emass import Functional

input = read_input("input.py", namespace = {"Escan": Functional})

G = array([0,0,0], dtype="float64")
X = array([0,0,1], dtype="float64")
structure = fill_structure(input.vff.lattice.cell, input.vff.lattice)
input.escan.nbstates = len(structure.atoms) * 4 + 4
input.escan.tolerance = 1e-12

orig = input.escan( structure, direction=(0,0,1), outdir="results/emass/100", \
                      do_relax_kpoint=False, type="e" )
assert abs(orig.mass[0] - 0.4381) < 0.01 # Gamma conduction emass 

result = input.escan( structure, direction=((0,0,1), (1,1,1)), outdir="results/hmass", \
                      do_relax_kpoint=False, type="h", bandgap=orig.extract_bg )
assert abs(result.mass[0][0] - 0.2769) < 0.01 # Gamma conduction heavy hmass (100 direction)
assert abs(result.mass[0][2] - 0.2059) < 0.01 # Gamma conduction light hmass (100 direction)
assert abs(result.mass[1][0] - 0.6885) < 0.01 # Gamma conduction heavy hmass (111 direction)
assert abs(result.mass[1][2] - 0.1460) < 0.01 # Gamma conduction light hmass (111 direction)
