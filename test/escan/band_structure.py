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

from numpy import array
from pylada.escan import read_input, exec_input, ReducedBPoints
from pylada.crystal.binary import zinc_blende

X = array( [1,0,0], dtype="float64" )
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" )
W = array( [0, 0.5,1], dtype="float64" )

input = read_input('input.py')
kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

structure = zinc_blende().to_structure(subs={'A':'Si', 'B':'Si'}) 
structure.scale = 5.45

kescan.fft_mesh = 14, 14, 14
kescan.kpoints = ReducedBPoints(density=20) + (X, G) + (G, L)
result = kescan( structure, outdir='results/kescan', 
                 nbstates = len(structure.atoms) * 4 + 4,
                 eref = None )
