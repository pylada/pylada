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

def functional():
  from pylada.vff.vff import Vff
  vff = Vff()
  vff["Si", "Si"] = 2.35905923, 45.375785351, -150.363195485, 519.230157133
  vff["Ge", "Ge"] = 2.44167023, 37.816422666, -124.830189020, 250.179861133
  vff["Si", "Ge"] = 2.39642606, 41.875535703, -145.994430120, 366.853558523
  vff["Si", "Si", "Si"] = "tet", 0, 12.579328566, 44.930324684, 359.422663897
  vff["Ge", "Ge", "Ge"] = "tet", 0, 10.383093493, 55.677487481
  vff["Si", "Ge", "Si"] = "tet", 0, 11.902809727, 20.193404352
  vff["Ge", "Si", "Ge"] = "tet", 0, 11.902809727, 20.193404352
  vff["Si", "Si", "Ge"] = "tet", 0, 12.245971457, 32.561864518, 179.706525419
  vff["Si", "Ge", "Ge"] = "tet", 0, 11.147853920, 37.930591733

  return vff

def crystal():
  from pylada.crystal import Structure

  structure = Structure( 10.0, 0.5, 0.5, 
                         0.00, 0.0, 0.5,
                         0.00, 0.5, 0.0, scale=5.5 )\
                       .add_atom(pos=(0.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(0.25, 0.25, 0.25), type="Ge") \
                       .add_atom(pos=(1.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(1.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(2.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(2.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(3.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(3.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(4.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(4.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(5.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(5.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(6.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(6.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(7.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(7.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(8.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(8.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(9.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(9.25, 0.25, 0.25), type="Si")  
  return structure

def test_relaxall(cartesian):
  """ Tests gradients as implemented in the functional. """
  from numpy import abs
  from pylada.vff.functional import Functional

  vff = functional()
  vff.cartesian = cartesian
  structure = crystal()
  epsilon = 1e-5
  
  # gets functions for current structure
  x, jacobian, energy, make_structure                                          \
       = Functional(copy=vff)._getfuncs_relaxall(structure)
  stress, forces = vff.jacobian(structure)

  jac0 = jacobian(x)
  for i in xrange(6):
    xminus = x.copy()
    xminus[i] = xminus[i] - epsilon * 0.5
    eminus = energy(xminus)
    xplus = x.copy()
    xplus[i] = xplus[i] + epsilon * 0.5
    eplus = energy(xplus)

    assert abs(jac0[i] - (eplus - eminus) / epsilon ) < 1e-6

  epsilon = 1e-7
  for i in xrange(6, len(x)):
    xminus = x.copy()
    xminus[i] = xminus[i] - epsilon * 0.5
    eminus = energy(xminus)
    xplus = x.copy()
    xplus[i] = xplus[i] + epsilon * 0.5
    eplus = energy(xplus)

    assert abs(jac0[i] - (eplus - eminus) / epsilon ) < 1e-6

def test_epi(direction, cartesian):
  from numpy import abs
  from pylada.vff.functional import Functional

  vff = functional()
  vff.cartesian = cartesian
  vff.direction = direction
  structure = crystal()
  
  # gets functions for current structure
  x, jacobian, energy, make_structure                                          \
       = Functional(copy=vff)._getfuncs_epi(structure)
  stress, forces = vff.jacobian(structure)

  jac0 = jacobian(x)
  epsilon = 1e-7
  for i in xrange(0, len(x)):
    xminus = x.copy()
    xminus[i] = xminus[i] - epsilon * 0.5
    eminus = energy(xminus)
    xplus = x.copy()
    xplus[i] = xplus[i] + epsilon * 0.5
    eplus = energy(xplus)

    assert abs(jac0[i] - (eplus - eminus) / epsilon ) < 1e-6

if __name__ == '__main__':
  test_relaxall(False)
  test_relaxall(True)
  test_epi([1, 0, 0], False)
  test_epi([1, 0, 0], True)
  test_epi([1, 2, 0], False)
  test_epi([1, 2, 0], True)
  test_epi([0, 2, 3], False)
  test_epi([0, 2, 3], True)
