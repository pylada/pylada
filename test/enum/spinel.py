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

from pylada.crystal.A2BX4 import b5
from pylada.enumeration import Enum
from spinel_check import checkme


# create spinel lattice and make it inverse.
lattice = Enum(b5())
for site in lattice.sites:
  if "X" in site.type: continue
  site.type = 'A' if 'B' in site.type else ('A', 'B')
lattice.find_space_group()

# look at inverse only.
def check_concentration(x, flavorbase, smith):
  from pylada.enumeration import as_numpy

  types = as_numpy(x, flavorbase) 
  result = 2*len(types[types == 1]) == len(types)

  return result, None if result else False

result = [i for i, dummy, dummy, dummy in lattice.xn(2)]
result.append([i for i, dummy, dummy, dummy in lattice.xn(4)])

# creates file with result.
# with open("spinel_check.py", "w") as file: 
#   file.write('checkme = {0}\n'.format(repr(result)))

for a, b in zip(checkme, result): assert a == b
