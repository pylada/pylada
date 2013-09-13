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

def test_outersum():
  from numpy import array, outer, all, abs, zeros
  from pylada.ce.cppwrappers import outer_sum
 
  for type in ['int8', 'int16', 'float', 'float64']:
    tensor = array([0], dtype=type)
    spins = array([3], dtype=type)
    outer_sum(tensor, spins) # slow way to do a sum...
    assert all(abs(spins-tensor) < 1e-8)
    outer_sum(tensor, spins) # slow way to do a sum...
    assert all(abs(2*spins-tensor) < 1e-8)

    tensor = array([0, 0], dtype=type)
    spins = array([3, 4], dtype=type)
    outer_sum(tensor, spins) # slow way to do a sum...
    assert all(abs(spins-tensor) < 1e-8)
    outer_sum(tensor, spins) # slow way to do a sum...
    assert all(abs(2*spins-tensor) < 1e-8)

    tensor = array([[0]], dtype=type)
    spins = array([3], dtype=type), array([-2], dtype=type)
    outer_sum(tensor, *spins)
    assert all(abs(outer(*spins) - tensor) < 1e-8)
    outer_sum(tensor, *spins)
    assert all(abs(2*outer(*spins) - tensor) < 1e-8)

    tensor = array([[0, 0]], dtype=type)
    spins = array([3], dtype=type), array([-2, 6], dtype=type)
    outer_sum(tensor, *spins)
    assert all(abs(outer(*spins) - tensor) < 1e-8)
    outer_sum(tensor, *spins)
    assert all(abs(2*outer(*spins) - tensor) < 1e-8)
    
    tensor = zeros((2,3), dtype=type)
    spins = array([3, 2], dtype=type), array([-2, 6, 4], dtype=type)
    outer_sum(tensor, *spins)
    assert all(abs(outer(*spins) - tensor) < 1e-8)
    outer_sum(tensor, *spins)
    assert all(abs(2*outer(*spins) - tensor) < 1e-8)

    tensor = zeros((2, 2, 3), dtype=type)
    spins = array([3, 2], dtype=type), array([10, 5], dtype=type), array([-2, 6, 4], dtype=type)
    outer_sum(tensor, *spins)
    assert all(abs(outer(outer(*spins[:2]), spins[-1]).reshape(2, 2, 3) - tensor) < 1e-8)
    outer_sum(tensor, *spins)
    assert all(abs(2*outer(outer(*spins[:2]), spins[-1]).reshape(2, 2, 3) - tensor) < 1e-8)

  # no args
  result = tensor.copy()
  outer_sum(tensor)
  assert all(abs(result - tensor) < 1e-8) # should do nothing

def test_product():
  from pylada.ce.cppwrappers import ProductILJ
  from pylada.crystal import Structure
  a = [u for u in ProductILJ(xrange(3), 1)]
  assert len(a) == 3
  assert a[0] == (0,)
  assert a[1] == (1,)
  assert a[2] == (2,)

  a = [u for u in ProductILJ(xrange(3), 2)]
  assert len(a) == 3
  assert a[0] == (0, 1)
  assert a[1] == (0, 2)
  assert a[2] == (1, 2)

  a = [u for u in ProductILJ(xrange(3), 3)]
  assert len(a) == 1
  assert a[0] == (0, 1, 2)

  a = [u for u in ProductILJ(xrange(3), 4)]
  assert len(a) == 0

  b = ('a', 5, 6j, Structure())
  a = [tuple([id(v) for v in u]) for u in ProductILJ(b, 2)]
  assert len(a) == 6
  assert a[0] == (id(b[0]), id(b[1]))
  assert a[1] == (id(b[0]), id(b[2]))
  assert a[2] == (id(b[0]), id(b[3]))
  assert a[3] == (id(b[1]), id(b[2]))
  assert a[4] == (id(b[1]), id(b[3]))
  assert a[5] == (id(b[2]), id(b[3]))

if __name__ == '__main__':
# test_outersum()
  test_product()
