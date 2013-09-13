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

from numpy import diag, array, dot
from numpy.linalg import norm, inv
from matplotlib.pyplot import scatter, show
from numpy import array
from operator import itemgetter
from pylada.crystal.binary import zinc_blende
from pylada.crystal import dnc_iterator

cell = array([[10.0, 0.5, 0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0]], dtype="float64")
cell = array([[-13. ,   0. ,   0.5], [ 13. ,   0. ,   0.5], [  0. ,  19. ,   0. ]])

structure = zinc_blende().to_structure(cell) #diag([10,1,1]))
# print structure
# print dnc_iterator(structure, 45, 1.5).mesh
result = [False for u in structure.atoms]
for box in dnc_iterator(structure, 15, 0.8):

  all = [u for u in box]
  small_box = [(j, structure.atoms[i].pos + trans) for j, (i, trans, f) in enumerate(all) if f == True]
  large_box = array([structure.atoms[i].pos + trans for i, trans, f in all])
  for index, trans, in_small_box in all:
    assert result[index] == False
    result[index] = False

  for index, pos in small_box:
    sizes = sorted([ (i, norm(u)) for i, u in enumerate(large_box - pos) ], key=itemgetter(1))
#   print structure.atoms[all[index][0]].pos, all[index][0]
#   for i, bondlength in sizes:
#     print " ", bondlength, structure.atoms[all[i][0]].pos, structure.atoms[all[i][0]].pos + all[i][1]
    for i in range(2, 5): 
      assert abs(sizes[i][1] - sizes[1][1]) < 1e-12
    assert abs(sizes[5][1] - sizes[1][1]) > 1e-12
    
    
