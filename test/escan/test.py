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

from sys import exit
from shutil import rmtree
from os.path import exists, join
from math import ceil, sqrt
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world
from pylada.vff import Vff
from pylada.crystal import Structure, Lattice, fill_structure
from pylada.escan import read_input

input = read_input("input.py")

structure = Structure()
structure.set_cell = (4, 0, 0.5),\
                     (0, 1,   0),\
                     (0, 0, 0.5)
structure = fill_structure(structure.cell)
for i, atom in enumerate(structure.atoms):
  atom.type = "Si" if i < len(structure.atoms)/2 else "Ge"


result_str = Structure()
result_str.scale = 5.450000e+00
result_str.set_cell = (4.068890e+00, -4.235770e-18, 5.083297e-01),\
                     (-1.694308e-17, 1.016103e+00, 2.238072e-18),\
                     (-2.252168e-03, 8.711913e-18, 5.083297e-01)
result_str.weight = 1.000000e+00
result_str.name = ""
result_str.energy = 0.0938967086716
result_str.add_atom = (0.000000e+00, 0.000000e+00, 0.000000e+00), "Si", 0,  0
result_str.add_atom = (2.541649e-01, 2.473273e-01, 2.541649e-01), "Si", 1,  0
result_str.add_atom = (3.567265e+00, 5.062000e-01, -8.956567e-03), "Si", 0,  0
result_str.add_atom = (3.821430e+00, 7.572301e-01, 2.452083e-01), "Si", 1,  0
result_str.add_atom = (3.065136e+00, -1.851371e-03, -1.515736e-02), "Si", 0,  0
result_str.add_atom = (3.319301e+00, 2.491787e-01, 2.390075e-01), "Si", 1,  0
result_str.add_atom = (2.563510e+00, 5.080514e-01, -2.186176e-02), "Si", 0,  0
result_str.add_atom = (2.817675e+00, 7.553787e-01, 2.323031e-01), "Si", 1,  0
result_str.add_atom = (2.055673e+00, -6.642716e-03, -2.235452e-02), "Ge", 0,  0
result_str.add_atom = (2.309838e+00, 2.539701e-01, 2.318104e-01), "Ge", 1,  0
result_str.add_atom = (1.539450e+00, 5.026981e-01, -1.446032e-02), "Ge", 0,  0
result_str.add_atom = (1.793614e+00, 7.607320e-01, 2.397046e-01), "Ge", 1,  0
result_str.add_atom = (1.024061e+00, -5.353269e-03, -7.401445e-03), "Ge", 0,  0
result_str.add_atom = (1.278226e+00, 2.526806e-01, 2.467634e-01), "Ge", 1,  0
result_str.add_atom = (5.078370e-01, 5.014086e-01, 4.927555e-04), "Ge", 0,  0
result_str.add_atom = (7.620018e-01, 7.620214e-01, 2.546576e-01), "Ge", 1,  0


X = array( [0,0,1], dtype="float64" ) # cartesian 2pi/a
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" ) 
W0 = array( [1, 0.5,0], dtype="float64" ) 
W1 = array( [1, 0,0.5], dtype="float64" ) 
W2 = array( [0, 1,0.5], dtype="float64" ) 

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
jobs = [\
         # at gamma, code uses Krammer degeneracy
         (G,   "VBM", -0.4, array([-0.47992312, -0.47992312,-0.67148097,-0.67148097])), 
         (G, "Gamma",  0.4, array([ 0.47368306,  0.47368306, 0.49199994, 0.49199994])), 
         (X,     "X",  0.4, array([ 0.51468608,  0.51479076, 0.5148467 , 0.5149207 ])),
         (L,     "L",  0.4, array([ 0.72789198,  0.72789198, 0.73165765, 0.73165765])),
         (W1,   "W1",  0.4, array([ 0.89170814,  0.89170822, 0.96097565, 0.96097601])),
         (W2,   "W2",  0.4, array([ 0.89174454,  0.89174462, 0.9608853 , 0.96088566]))
       ]

# end of input 


# computes vff results and checks
# if exists("work") and world.rank == 0: rmtree("work")
vff_out = input.vff(structure, outdir = "results", comm = world)
exit(0)
# if world.rank == 0:
#   solo = vff_out.solo()
#   diff = 0e0
#   for a, b in zip(solo.structure.atoms, result_str.atoms):
#     assert a.type == b.type
#     assert a.site == b.site
#     diff += dot(a.pos-b.pos, a.pos-b.pos)

#   assert diff / float(len(structure.atoms)) < 1e-8, diff 
#   assert abs(solo.energy - result_str.energy) < 1e-8, abs(solo.energy - result_str.energy) 

# some kpoints
# launch escan for different jobs.
for kpoint, name, ref, expected_eigs in jobs:
  out = input.escan( structure, join("results", name), comm = world,\
               kpoint = kpoint,
               eref = None,\
               nbstates = 26,\
               workdir="work")
  if world.rank == 0:
    solo = out.solo()
    print solo.success, solo.nnodes
    print solo.eigenvalues
    print solo.convergence
  break;
  # if world.rank == 0: print "Ok - %s: %s -> %s: %s" % (name, kpoint, escan.kpoint, eigenvalues)
  
