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

def test_checkkeywords():
  """ Checks that check keyword works. """
  from pylada.ce._single_site_factory import check_keywords
  from pylada.error import KeyError, ValueError, TypeError

  check_keywords(True, False)
  check_keywords(False, True)
  check_keywords(False, False, B2=2)
  check_keywords(True, False, B2=3, B8=5, B15=12)

  try: check_keywords(False, False)
  except ValueError: pass
  else: raise Exception()

  try: check_keywords(False, False, B2=0)
  except ValueError: pass
  else: raise Exception()
  
  try: check_keywords(False, False, B2=-1)
  except ValueError: pass
  else: raise Exception()
  
  try: check_keywords(True, False, b2=1)
  except KeyError: pass
  else: raise Exception()

  try: check_keywords(True, False, B1=1)
  except KeyError: pass
  else: raise Exception()

  try: check_keywords(True, False, B0=1)
  except KeyError: pass
  else: raise Exception()

  try: check_keywords(True, False, B5='a')
  except TypeError: pass
  else: raise Exception()

def test_pairs():
  """ Tests pair generation. """
  from numpy import all, sum
  from pylada.ce._single_site_factory import factory
  from pylada.crystal import bravais
  from pylada.error import ValueError

  lattice = bravais.fcc()
  try: results = factory(lattice, B2=1)
  except ValueError: pass
  else: raise Exception()

  lattice[0].type = ['A', 'B']
  results = factory(lattice, B2=1)
  assert len(results) == 1
  assert all(results[0].spins['position'][0] == [0,0,0])
  assert abs(sum(results[0].spins['position'][1]**2) - 2*0.5*0.5) < 1e-8
  assert len(results[0]._symmetrized) == 12

  results = factory(lattice, B2=2)
  assert len(results) == 2
  assert all(results[0].spins['position'][0] == [0,0,0])
  assert sum(results[0].spins['position'][1]**2) == 2*0.5*0.5
  assert len(results[0]._symmetrized) == 12
  assert all(results[1].spins['position'][0] == [0,0,0])
  assert abs(sum(results[1].spins['position'][1]**2) - 1e0) < 1e-8
  assert len(results[1]._symmetrized) == 6

  results = factory(lattice, B2=3)
  assert len(results) == 3
  assert all(results[0].spins['position'][0] == [0,0,0])
  assert sum(results[0].spins['position'][1]**2) == 2*0.5*0.5
  assert len(results[0]._symmetrized) == 12
  assert all(results[1].spins['position'][0] == [0,0,0])
  assert abs(sum(results[1].spins['position'][1]**2) - 1e0) < 1e-8
  assert len(results[1]._symmetrized) == 6
  assert all(results[2].spins['position'][0] == [0,0,0])
  assert abs(sum(results[2].spins['position'][1]**2) - 1-2*0.25) < 1e-8
  assert len(results[2]._symmetrized) == 24

def test_3b():
  """ Tests triplet generation. """
  from pylada.ce._single_site_factory import factory
  from pylada.crystal import bravais

  lattice = bravais.fcc()
  lattice[0].type = ['A', 'B']
  results = factory(lattice, B3=1)
  assert len(results) == 4
  results = factory(lattice, B3=2)
  assert len(results) == 9

if __name__ == '__main__':
  test_checkkeywords()
  test_pairs()
  test_3b()
