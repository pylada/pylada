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

from pickle import loads, dumps
def test_molecule():
  from numpy import all, abs, array
  from pylada.dftcrystal.molecule import Molecule

  l = ['Ge', 'Si']
  a = Molecule()
  a.add_atom(0, 0, 0, 'Si')
  a.add_atom(1, 2, 3, l)

  c, d = loads(dumps((l, a)))
  assert len(d.atoms) == 2
  assert d.__class__ is Molecule
  assert d.keyword == a.keyword
  assert all(abs(d.atoms[0].pos - a.atoms[0].pos) < 1e-8)
  assert all(abs(d.atoms[1].pos - a.atoms[1].pos) < 1e-8)
  assert all(d.atoms[0].type == a.atoms[0].type)
  assert all(d.atoms[1].type is c)
  assert all(array(c) == l)

def test_hamiltonian():
  from pylada.dftcrystal.hamiltonian import Dft

  a = Dft()
  b = loads(dumps(a))
  assert repr(a) == repr(b)

  a.b3lyp = True
  a.xxlgrid = True
  assert repr(a) == repr(loads(dumps(a)))

def test_electronic():
  from pylada.dftcrystal.electronic import Electronic

  a = Electronic()
  assert repr(a) == repr(loads(dumps(a)))
  a.shrink = [8, 8, 9], None
  a.tolinteg = 8, 8, 8, 8, 14
  a.dft.b3lyp = True
  a.dft.exchange = 'lda'
  a.dft.xxlgrid = True
  assert repr(a) == repr(loads(dumps(a)))

def test_optgeom():
  from pylada.dftcrystal.optgeom import OptGeom
  from pylada.dftcrystal.input import AttrBlock
  a = OptGeom()
  assert repr(a) == repr(loads(dumps(a)))
  a.maxcycle = 50
  a.fulloptg = True
  assert repr(a) == repr(loads(dumps(a)))
  assert repr(OptGeom()) != repr(loads(dumps(a)))

  c = AttrBlock()
  c.optgeom = a
  c.optgeom.maxcycle = 10
  assert repr(c) == repr(loads(dumps(c)))

def test_basis():
  from pylada.dftcrystal.basis import BasisSet, Shell

  a = BasisSet()
  assert repr(a) == repr(loads(dumps(a)))
  a['H'] = [ Shell( 's', 1.0, 
                    a0=[18.731137, 0.0334946],
                    a1=[2.8253937, 0.23472695],
                    a2=[0.6401217, 0.81375733] ),
             Shell('s', 0.0, a0=[0.1612778, 1.0]), 
             Shell('p', 0.0, a0=[1.1, 1.0]) ]
  assert repr(a) == repr(loads(dumps(a)))
  assert repr(BasisSet()) != repr(loads(dumps(a)))

def test_functional():
  from pylada.dftcrystal import Functional, Shell
  from pylada.dftcrystal.input import print_input
  
  a = Functional()
  assert repr(a) == repr(loads(dumps(a)))
  a.basis['H'] = [ Shell( 's', 1.0, 
                          a0=[18.731137, 0.0334946],
                          a1=[2.8253937, 0.23472695],
                          a2=[0.6401217, 0.81375733] ),
                   Shell('s', 0.0, a0=[0.1612778, 1.0]), 
                   Shell('p', 0.0, a0=[1.1, 1.0]) ]
  a.shrink = [8, 8, 9], None
  a.tolinteg = 8, 8, 8, 8, 14
  a.dft.b3lyp = True
  a.dft.exchange = 'lda'
  a.dft.xxlgrid = True
  a.optgeom.maxcycle = 50
  a.optgeom.fulloptg = True
  assert repr(Functional()) != repr(loads(dumps(a)))
  assert repr(a) == repr(loads(dumps(a)))
  b = loads(dumps(a))
  assert print_input(a.output_map(crystal=a)) == print_input(b.output_map(crystal=b))

if __name__ == '__main__':
  test_molecule()
  test_hamiltonian()
  test_electronic()
  test_optgeom()
  test_basis()
  test_functional()
