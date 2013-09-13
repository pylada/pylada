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
  import gc
  from numpy import all, abs, ones
  from sys import getrefcount
  from pylada.crystal.binary import zinc_blende
  from pylada.error import ValueError, AttributeError
  from pylada.vff import Node

  structure = zinc_blende()
  assert getrefcount(structure[0]) == 2 # structure, getrefcount arg.
  nodeA = Node(structure[0], 0)
# assert getrefcount(structure[0]) == 3 # structure, node, getrefcount arg.

# assert nodeA.center is structure[0]
# assert all(abs(nodeA.pos - structure[0].pos) < 1e-8)
# assert nodeA.type == structure[0].type
# nodeA.pos[0] = 0.5
# nodeA.type = 'C'
# assert all(abs(nodeA.pos - structure[0].pos) < 1e-8)
# assert nodeA.type == structure[0].type
# assert getrefcount(structure[0]) == 3 # structure, node, getrefcount arg.

# assert len(nodeA) == 0
# assert getrefcount(structure[1]) == 2 # structure, getrefcount arg.
# nodeB = Node(structure[1], 0)
# assert len(nodeB) == 0
# assert getrefcount(structure[0]) == 3 # structure, node, getrefcount arg.
# assert getrefcount(structure[1]) == 3 # structure, node, getrefcount arg.

# # should increase ref count on node by 1
# assert getrefcount(nodeA) == 2
# assert getrefcount(nodeB) == 2
# assert nodeA.link(nodeB, [0,0,0])
# assert len(nodeA) == 1
# assert len(nodeB) == 1
# assert nodeA[0][0] is nodeB
# assert nodeB[0][0] is nodeA
# assert all(abs(nodeA[0][1] + nodeB[0][1]) < 1e-8)
# assert all(abs(nodeA[0][1]) < 1e-8)
# # should not increase ref count on atom
# assert getrefcount(structure[0]) == 3 # structure, node, getrefcount arg.
# # should increase ref count on node by 1
# assert getrefcount(nodeA) == 3
# assert getrefcount(nodeB) == 3

# assert not nodeA.link(nodeB, [0,0,0])
# assert len(nodeB) == 1 and len(nodeA) == 1
# assert not nodeA.link(nodeB)
# assert len(nodeB) == 1 and len(nodeA) == 1
# assert not nodeB.link(nodeA, [0,0,0])
# assert len(nodeB) == 1 and len(nodeA) == 1
# assert not nodeB.link(nodeA)
# assert len(nodeB) == 1 and len(nodeA) == 1
# assert getrefcount(nodeA) == 3
# assert getrefcount(nodeB) == 3

# try: nodeA.link(nodeA)
# except ValueError: pass
# else: raise Exception()
# assert getrefcount(nodeA) == 3
# assert getrefcount(nodeB) == 3

# # ref count should increase only be 1, despite both endpoints pointing to
# # nodeA
# assert nodeA.link(nodeA, [1, 0, 0])
# assert len(nodeA) == 2
# assert len(nodeB) == 1
# assert getrefcount(nodeA) == 4 
# assert getrefcount(nodeB) == 3

# # check iteration
# iter = nodeA.__iter__()
# a, b = iter.next()
# assert a is nodeB
# assert all(abs(b) < 1e-8)
# a, b = iter.next()
# assert a is nodeA
# assert all(abs(b-[1, 0, 0]) < 1e-8)
# try: iter.next()
# except StopIteration: pass
# else: raise Exception()
# del iter
# del a
# del b

# # check deletion
# nodeB.clear()
# gc.collect()
# assert len(nodeA) == 1
# assert len(nodeB) == 0
# assert getrefcount(nodeA) == 3 
# assert getrefcount(nodeB) == 2
# # check deletion
# nodeA.clear()
# gc.collect()
# assert getrefcount(nodeA) == 2
# assert len(nodeA) == 0

# # check trnaslation is opposite for opposite endpoints.
# nodeB.link(nodeA, [1, 0, 0])
# assert len(nodeA) == 1
# assert len(nodeB) == 1
# assert all(abs(nodeB[0][1] - [1, 0, 0]) < 1e-8)
# assert all(abs(nodeA[0][1] + nodeB[0][1]) < 1e-8)
# assert not nodeA.link(nodeB, [-1, 0, 0])
# 
# # check we can set gradient
# try: nodeA.gradient
# except AttributeError: pass
# else: raise Exception()
# nodeA.gradient = ones(3, dtype='float64')
# assert all(abs(nodeA.gradient-1e0) < 1e-8)
# del nodeA.gradient
# try: nodeA.gradient
# except AttributeError: pass
# else: raise Exception()


if __name__ == '__main__':
  test()
