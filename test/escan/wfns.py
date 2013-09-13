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

""" Checks that we know how to go from g-space to r-space and back.

    Uses two methods: (i) *manual* Fourier transform using numpy tensorial
    products, (ii) d3fft from escan fortran code. The results from these two
    methods are compare to one another.
""" 
import pickle
from sys import exit
from os import getcwd
from os.path import join
from numpy import matrix, array
from numpy.linalg import norm
from boost.mpi import world
from pylada.opt import read_input, redirect
from pylada.opt.changedir import Changedir
from pylada.escan import Escan, soH, bandgap
from pylada.escan._escan import to_realspace
from pylada.escan._wfns import rtog_fourrier
from pylada.vff import Vff
from pylada.crystal import fill_structure, sort_layers, FreezeCell, nb_valence_states

# reads input file.
global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
input = read_input("input.py", global_dict=global_dict)

# creating unrelaxed structure.
input.vff.lattice.set_as_crystal_lattice()
# cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
# structure = sort_layers(fill_structure(cell), array([0,0,2.5]))
# for i, atom in enumerate(structure.atoms):
#   atom.type = "Si" if i % 10 < 6 else "Ge"
structure = fill_structure(input.escan.vff.lattice.cell)
structure.atoms[0].type = "Si"
structure.scale = 5.55
structure.freeze = FreezeCell.a0 | FreezeCell.a1
input.escan.fft_mesh  = 14, 14, 14

out = input.escan( structure,\
                   outdir=join("results", "wfns"),\
                   comm=world )

with Changedir(out.directory) as directory:
  nbstates = out.escan.nbstates if out.escan.potential == soH and norm(out.escan.kpoint)\
             else out.escan.nbstates / 2
  out.raw_gwfns
  out.escan._write_incar(out.comm, out.structure)
  with redirect(fout="") as streams:
    rwfns, rvectors = to_realspace(input.escan, out.raw_gwfns, world)
assert out.rvectors.shape == rvectors.shape
for a, b in zip(out.rvectors, rvectors):
  assert norm(a-b) < 1e-8

assert out.raw_rwfns[0].shape == rwfns.shape
for a, b in zip(out.raw_rwfns[0].flat, rwfns.flat):
  assert abs(a-b) < 1e-8

gwfns = rtog_fourrier(rwfns, rvectors, out.gvectors, out.comm)
assert gwfns.shape == out.raw_gwfns.shape
for a, b in zip(out.raw_gwfns.flat, gwfns.flat):
  assert abs(a-b) < 1e-8
