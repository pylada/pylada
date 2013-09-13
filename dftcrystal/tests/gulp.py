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
  """ Tests writing out gulp files. """
  from numpy import array, all, abs
  from pylada.dftcrystal import Crystal
  from pylada.crystal.write import gulp

  a = Crystal(136, 4.63909875, 2.97938395, \
              ifhr=0, \
              shift=0)\
             .add_atom(0, 0, 0, 'Ti')\
             .add_atom(0.306153, 0.306153, 0, 'O')

  string = gulp(a).splitlines()
  assert string[0] == 'vectors'
  assert all(abs(array(string[1].split(), dtype='float64') - [4.63909875, 0, 0]) < 1e-5) 
  assert all(abs(array(string[2].split(), dtype='float64') - [0, 4.63909875, 0]) < 1e-5) 
  assert all(abs(array(string[3].split(), dtype='float64') - [0, 0, 2.97938395]) < 1e-5) 
  assert string[4] == 'spacegroup'
  assert string[5] == '136'
  assert string[6].rstrip().lstrip() == 'cartesian'
  assert string[7].split()[:2] == ['Ti', 'core']
  assert all(abs(array(string[7].split()[2:], dtype='float64')) < 1e-5)
  assert string[8].split()[:2] == ['O', 'core']
  assert all(abs(array(string[8].split()[2:], dtype='float64') - [1.420274, 1.420274, 0]) < 1e-5)

if __name__ == '__main__': 
  test()





