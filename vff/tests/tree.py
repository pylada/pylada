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

def test_tree():
  from numpy import all, array, dot, sum, any
  from pylada.crystal import binary, supercell
  from pylada.vff import build_tree
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = build_tree(a, overlap=0.5)
  
  for center in b:
    positions = []
    for i, (bond, vector) in enumerate(center):
      position = bond.pos + dot(a.cell, vector)
      assert abs(sum((position - center.pos)**2) - 0.25*0.25*3) < 1e-8
      assert all( [any(abs(array(p) - position[None, :]) > 1e-8) for p in positions] )
      positions.append(position)
    assert i == 3

def test_disorder(lim=8):
  from numpy import all, array, dot
  from numpy.random import random, randint
  from numpy.linalg import det
  from pylada.crystal import binary, supercell
  from pylada.vff import build_tree


  lattice = binary.zinc_blende()
  for i in xrange(10):
    cell = randint(-lim, lim, (3,3))
    while det(cell) == 0: cell = randint(-lim, lim, (3,3))
    a = supercell(lattice, dot(lattice.cell, cell))
  
    b = build_tree(a, overlap=0.8)
    ids = [id(node.center) for node in b]
    connections = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
  
    epsilon = random((3,3)) * 0.1
    epsilon = epsilon + epsilon.T
    a.cell += dot(epsilon, a.cell)
    for atom in a: atom.pos += dot(epsilon, atom.pos)
    
    b = build_tree(a, overlap=0.8)
    c = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
    assert all(connections == c)
    
    b = build_tree(a, overlap=0.8)
    for atom in a: atom.pos += random(3) * 0.05 - 0.025
    c = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
    assert all(connections == c)
  return a

if __name__ == '__main__':
  test_tree()
  structure = test_disorder()
