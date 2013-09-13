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

from pylada.crystal import Structure
from pylada.pcm import Clj, bond_name
from pylada.physics import a0, Ry
from quantities import angstrom, eV, hartree

clj  = Clj()
""" Point charge + r^12 + r^6 model. """
clj.ewald_cutoff = 80 * Ry

clj.charges["A"] = -1.0
clj.charges["B"] =  1.0

structure = Structure()
structure.set_cell = (1,0,0),\
                     (0,1,0),\
                     (0,0,1)
structure.scale = 50
structure.add_atom = (0,0,0), "A"
structure.add_atom = (a0.rescale(angstrom)/structure.scale,0,0), "B"

print clj.ewald(structure).energy, hartree.rescale(eV)


from pylada.crystal.A2BX4 import b5
from pylada.crystal import fill_structure
from numpy import array
clj.ewald_cutoff = 20 * Ry
lattice = b5()
lattice.sites[4].type='A'
structure = fill_structure(lattice.cell, lattice)
structure.scale = 8.0

clj.charges["A"] =  3.0
clj.charges["B"] =  2.0
clj.charges["X"] = -2.0
print clj.ewald(structure).energy, -498.586
