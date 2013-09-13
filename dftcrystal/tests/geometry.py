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

def test_crystal():
  from numpy import array
  from pylada.dftcrystal.crystal import Crystal
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O')
  assert structure.print_input() is not None
  assert len(structure.print_input().split('\n')) == 9
  assert all(abs(array(structure.print_input().split()[1:-2], dtype='float64') \
             - [0, 0, 0, 136, 4.63909875, 2.97938395, 2, 22, 0, 0, 0, 8, 0.306153, 0.306153, 0]) < 1e-8)
  assert structure.print_input().split()[0] == 'CRYSTAL'
  assert structure.print_input().split()[-2] == 'END'
  assert structure.print_input().split()[-1] == 'CRYSTAL'

def test_supercell():
  from numpy import array, all, abs
  from pylada.dftcrystal import Supercell, Crystal
  from pylada.error import ValueError
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O') \
                  .append(Supercell([[2,0,0],[0,1,0],[0,0,1]]))
  assert len(structure.eval()) == 2 * 6
  structure[-1].matrix[2,2] = 2
  assert len(structure.eval()) == 2 * 2 * 6
  map = structure[-1].output_map()
  assert 'supercel' in map
  a = Supercell()
  a.read_input(map['supercel'])
  assert all(abs(a.matrix - structure[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Supercell': Supercell}).matrix - a.matrix) < 1e-8)

  # test errors
  try: structure[-1].matrix = array( [[0]*3]*3 )
  except ValueError: pass
  else: raise Exception()
  structure[-1].matrix = array([[1,0],[0,1]])
  try: structure.print_input()
  except ValueError: pass
  else: raise Exception()

  map = structure[-1].output_map()
  assert 'supercel' in map
  a = Supercell()
  a.read_input(map['supercel'])
  assert all(abs(a.matrix - structure[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Supercell': Supercell}).matrix - a.matrix) < 1e-8)


def test_elastic():
  """ Test elastic keyword. """
  from pickle import loads, dumps
  from numpy import array, all, abs, identity, dot
  from pylada.dftcrystal import Crystal, Elastic

  epsilon = array([[0, 0, 0], [0, 0, 0], [0, 0, 0.1]])
  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  structure = crystal.eval()
  crystal.append(Elastic([[0, 0, 0], [0, 0, 0], [0, 0, 0.1]]))
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  assert abs(crystal.eval()[0].pos[1] - crystal.eval()[0].pos[0]) < 1e-8
  assert abs(crystal.eval()[1].pos[2] / crystal.eval()[1].pos[1] - 1.1) < 1e-8
  assert abs(crystal.eval()[1].pos[1] - crystal.eval()[1].pos[0]) < 1e-8
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Elastic': Elastic}).matrix - a.matrix) < 1e-8)
  assert eval(repr(a), {'Elastic': Elastic}).is_epsilon == a.is_epsilon
  assert loads(dumps(a)).is_epsilon == a.is_epsilon
  assert all(abs(loads(dumps(a)).matrix - a.matrix) < 1e-8)

  z = dot(epsilon, structure.cell)
  crystal[-1].matrix = z
  crystal[-1].is_epsilon = False
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Elastic': Elastic}).matrix - a.matrix) < 1e-8)
  assert eval(repr(a), {'Elastic': Elastic}).is_epsilon == a.is_epsilon
  assert loads(dumps(a)).is_epsilon == a.is_epsilon
  assert all(abs(loads(dumps(a)).matrix - a.matrix) < 1e-8)

  epsilon = array([[0, 0, 0], [0, 0, 0], [0, 0.1, 0]])
  crystal[-1].matrix = epsilon
  crystal[-1].is_epsilon = True
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Elastic': Elastic}).matrix - a.matrix) < 1e-8)
  assert eval(repr(a), {'Elastic': Elastic}).is_epsilon == a.is_epsilon
  assert loads(dumps(a)).is_epsilon == a.is_epsilon
  assert all(abs(loads(dumps(a)).matrix - a.matrix) < 1e-8)

  z = dot(epsilon, structure.cell)
  crystal[-1].matrix = z
  crystal[-1].is_epsilon = False
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)
  assert all(abs(eval(repr(a), {'Elastic': Elastic}).matrix - a.matrix) < 1e-8)
  assert eval(repr(a), {'Elastic': Elastic}).is_epsilon == a.is_epsilon
  assert loads(dumps(a)).is_epsilon == a.is_epsilon
  assert all(abs(loads(dumps(a)).matrix - a.matrix) < 1e-8)



if __name__ == '__main__':
  test_crystal()
  test_supercell()
  test_elastic()
