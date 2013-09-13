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

def test_J0():
  from pylada.ce import cluster_factory
  from pylada.crystal import binary

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, J0=True)
  assert len(a) == 1
  assert a[0].order == 0

def test_J1():
  from numpy import all
  from pylada.ce import cluster_factory
  from pylada.crystal import binary

  # test multi-lattice with different occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, J1=True)
  assert len(a) == 2
  assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
  assert a[0].spins['sublattice'] == 0
  assert a[1].spins['sublattice'] == 1

  # test multi-lattice with same occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']

  a = cluster_factory(lattice, J1=True)
  assert len(a) == 1
  assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
  assert a[0].spins['sublattice'] == 0

def test_B2():
  from numpy import all
  from pylada.ce import cluster_factory
  from pylada.crystal import binary

  # test multi-lattice with different occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, B2=1)
  assert len(a) == 2
  assert all(a[0].spins['sublattice'] == [0, 1])
  assert all(abs(a[0].spins['position'][0]) < 1e-8)
  vector = a[0].spins[1]
  assert abs(sum((lattice[vector[1]].pos + vector[0])**2) - 3*0.25*0.25) < 1e-8
  assert all(a[1].spins['sublattice'] == [1, 0])
  assert all(abs(a[1].spins['position'][0]) < 1e-8)
  assert all(abs(a[1].spins['position'][1]) < 1e-8)
  a = cluster_factory(lattice, B2=2)
  assert len(a) == 4
  a = cluster_factory(lattice, B2=3)
  assert len(a) == 6

  # test multi-lattice with same occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']

  a = cluster_factory(lattice, B2=1)
  assert len(a) == 1
  a = cluster_factory(lattice, B2=2)
  assert len(a) == 2
  a = cluster_factory(lattice, B2=3)
  assert len(a) == 3

def test_B3():
  from numpy import all
  from pylada.ce import cluster_factory
  from pylada.crystal import binary

  def topos(s):
    return lattice[s[1]].pos + s[0]
  # test multi-lattice with different occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, B3=1)
  assert len(a) == 2
  for cluster in a:
    assert all(abs(abs(topos(cluster.spins[1]) - topos(cluster.spins[0])) - 0.25) < 1e-8)
    assert all(abs(abs(topos(cluster.spins[2]) - topos(cluster.spins[0])) - 0.25) < 1e-8)

  # test multi-lattice with same occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']
  a = cluster_factory(lattice, B3=1)
  assert len(a) == 1
  for cluster in a:
    assert all(abs(abs(topos(cluster.spins[1]) - topos(cluster.spins[0])) - 0.25) < 1e-8)
    assert all(abs(abs(topos(cluster.spins[2]) - topos(cluster.spins[0])) - 0.25) < 1e-8)

if __name__ == '__main__': 
  test_J0()
  test_J1()
  test_B2()
  test_B3()
