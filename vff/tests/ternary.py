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
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff

def test_inas():
  from numpy import identity, abs, all, dot, array
  from pylada.crystal.binary import zinc_blende
  from quantities import eV, angstrom

  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 6.5 #2.62332 * 2 / sqrt(3)  / 0.529177

  out = vff._pyeval(structure)
  assert abs(out.energy - 0.34958768908 * eV) < 1e-8
  assert abs(out.energy - vff.energy(structure)) < 1e-8
  assert all(abs(out.stress - vff.jacobian(structure)[0]) < 1e-8)
  assert all(abs(out.stress - identity(3) * -0.04096678 * eV/angstrom**3) < 1e-8)
  assert all(abs(out[0].gradient) < 1e-8)
  assert all(abs(out[1].gradient) < 1e-8)

  epsilon = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(epsilon, structure.cell)
  for atom in structure: atom.pos = dot(epsilon, atom.pos)
  out = vff._pyeval(structure)
  assert abs(out.energy - 0.527010806043 * eV) < 1e-8
  assert abs(out.energy - vff.energy(structure)) < 1e-8
  assert all(abs(out.stress - [[ -2.50890474e-02,  -2.95278697e-02,  0],
                               [ -2.95278697e-02,  -2.50890474e-02,  0],
                               [ 0, 0,  -1.85427515e-02]] * eV / angstrom**3) < 1e-6)
  assert all(abs(out.stress - vff.jacobian(structure)[0]) < 1e-8)
  assert all(abs(out[0].gradient - [0, 0, 1.09205526] * eV / angstrom) < 1e-6)
  assert all(abs(out[1].gradient - [0, 0, -1.09205526] * eV / angstrom) < 1e-6)
  assert all(abs([u.gradient for u in out] - vff.jacobian(structure)[1].magnitude) < 1e-8)

def test_ingaas():
  from numpy import abs, all, dot, array
  from quantities import eV, angstrom
  from pylada.crystal import Structure

  vff = functional()

  structure = Structure( 10.0, 0.5, 0.5, 
                         0.00, 0.0, 0.5,
                         0.00, 0.5, 0.0, scale=6.5 )\
                       .add_atom(pos=(0.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(0.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(1.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(1.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(2.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(2.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(3.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(3.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(4.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(4.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(5.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(5.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(6.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(6.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(7.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(7.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(8.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(8.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(9.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(9.25, 0.25, 0.25), type="As")  

  epsilon = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(epsilon, structure.cell)
  for atom in structure: atom.pos = dot(epsilon, atom.pos)
  out = vff._pyeval(structure)
  assert abs(out.energy - 12.7962141476*eV) < 1e-8
  assert abs(vff.energy(structure) - 12.7962141476*eV) < 1e-8
  stress = array([[ 0.07050804,  0.04862879,  0.00025269],
                  [ 0.04862879,  0.07050804, -0.00025269],
                  [ 0.00025269, -0.00025269,  0.06073765]]) * eV/angstrom**3
  assert all(abs(out.stress - stress) < 1e-6)
  assert all(abs(vff.jacobian(structure)[0] - stress) < 1e-6)
  gradients = array([u.gradient for u in out])
  check_gradients = array( [[  9.15933995e-16,  -5.96744876e-15,   1.12373933e+00],
                            [  2.22044605e-16,  -9.43689571e-15,  -1.12373933e+00],
                            [  6.19921859e-02,  -6.19921859e-02,   1.00514227e+00],
                            [ -5.99520433e-15,   2.19269047e-15,  -1.12373933e+00],
                            [  7.66544525e-02,  -7.66544525e-02,   1.02138259e+00],
                            [  3.82409256e-01,  -3.82409256e-01,  -1.65478066e+00],
                            [ -1.89756333e-02,   1.89756333e-02,   1.15847492e+00],
                            [ -4.16333634e-17,   8.81239526e-16,  -1.09205526e+00],
                            [ -7.54665596e-02,   7.54665596e-02,   1.14107325e+00],
                            [ -3.64621516e-01,   3.64621516e-01,  -5.74094841e-01],
                            [  7.66544525e-02,  -7.66544525e-02,   1.02138259e+00],
                            [  3.82409256e-01,  -3.82409256e-01,  -1.65478066e+00],
                            [ -1.89756333e-02,   1.89756333e-02,   1.15847492e+00],
                            [ -1.19973476e-14,  -2.32591724e-14,  -1.09205526e+00],
                            [ -1.37458746e-01,   1.37458746e-01,   1.25967031e+00],
                            [ -3.64621516e-01,   3.64621516e-01,  -5.74094841e-01],
                            [ -7.77156117e-16,  -3.10862447e-15,   1.12373933e+00],
                            [  8.29891711e-15,   0.00000000e+00,  -1.12373933e+00],
                            [  3.33066907e-15,   7.16093851e-15,   1.12373933e+00],
                            [ -4.44089210e-16,   7.85482790e-15,  -1.12373933e+00]])
  assert all(abs(gradients - check_gradients) < 1e-8)
  assert all(abs(vff.jacobian(structure)[1].magnitude - check_gradients) < 1e-8)

if __name__ == '__main__': 
  test_inas()
  test_ingaas()
