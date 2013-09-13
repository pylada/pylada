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

def test_single_counting():
  from pylada.crystal import binary, supercell
  from pylada.vff import build_tree
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[4, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = build_tree(a, overlap=0.5)
  
  n = 0
  for center in b:
    for endpoint, vector in center.sc_bond_iter():
      n += 1
      for other, v in endpoint.sc_bond_iter(): assert other is not center
      assert id(center) in [id(c) for c, v in endpoint]
  assert n == 2 * len(a)

def test_angle():
  from pylada.crystal import binary, supercell
  from pylada.vff import build_tree

  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[4, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = build_tree(a, overlap=0.5)

  for center in b:
    ids = [id(u.center) for u, d in center]
    angles = set([(id(u.center), id(v.center)) for (u, d), (v, d) in center.angle_iter()])
    for i, ida in enumerate(ids[:-1]):
      for idb in ids[i+1:]:
        if (ida, idb) in angles: assert (idb, ida) not in angles
        else: assert (idb, ida) in angles
    assert len(angles) == 6

if __name__ == '__main__':
  test_single_counting()
  test_angle()
