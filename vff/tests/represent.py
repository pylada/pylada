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

def test_rep_vff():
  from numpy import all, abs
  from pylada.vff.vff import Vff
  from pylada.vff import exec_input

  a = Vff()
  input = exec_input(repr(a))
  assert len(input.vff._parameters) == 0

  a = functional()
  input = exec_input(repr(a))
  assert ("In", "As") in input.vff
  assert all(abs(input.vff["In", "As"] - a["In", "As"]) < 1e-8)
  assert all(abs(input.vff["In", "As", 'Ga'] - a["Ga", "As", 'In']) < 1e-8)

def test_pickle_vff():
  from pickle import loads, dumps
  from numpy import all, abs
  from pylada.vff.vff import Vff

  a = Vff()
  b = loads(dumps(a)) 
  assert len(b._parameters) == 0

  a = functional()
  b = loads(dumps(a)) 
  assert ("In", "As") in b
  assert all(abs(b["In", "As"] - a["In", "As"]) < 1e-8)
  assert all(abs(b["In", "As", 'Ga'] - a["Ga", "As", 'In']) < 1e-8)

if __name__ == '__main__':
  test_rep_vff()
  test_pickle_vff()
