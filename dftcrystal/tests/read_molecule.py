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

def test():
  from numpy import abs, all
  from pylada.dftcrystal.parse import parse
  from pylada.dftcrystal import Molecule
  string = '\nMOLECULE\n1\n2\n1 0.0 -2.91352558499E-15 12.3754696347\n'        \
           '1 0.0 -2.91352558499E-15 13.1454696347\nEND MOLECULE'
  tree = parse(string)['']
  molecule = Molecule()
  molecule.read_input(tree['MOLECULE'])
  assert molecule.symmgroup == 1
  assert len(molecule) == 0
  assert len(molecule.atoms) == 2
  assert molecule.atoms[0].type == 'H'
  assert molecule.atoms[1].type == 'H'
  assert all(abs(molecule.atoms[0].pos - [0, 0, 12.3754696347]) < 1e-6)
  assert all(abs(molecule.atoms[1].pos - [0, 0, 13.1454696347]) < 1e-6)

if __name__ == '__main__': test()
