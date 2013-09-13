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
  from numpy import all, abs, array
  from pylada.dftcrystal.crystal import Crystal
  a = Crystal()
  a.raw = '0 0 0\n'\
          '136\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.symmgroup == 136
  assert a.shift == 0
  assert a.ifhr == 0
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'Ti' and all(abs(a.atoms[0].pos) < 1e-8)
  assert a.atoms[1].type == 'O' and all(abs(a.atoms[1].pos-[0.30615265, 0.30615265, 0.0]) < 1e-8)
  a.raw = a.raw
  assert a.symmgroup == 136
  assert a.shift == 0
  assert a.ifhr == 0
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'Ti' and all(abs(a.atoms[0].pos) < 1e-8)
  assert a.atoms[1].type == 'O' and all(abs(a.atoms[1].pos-[0.30615265, 0.30615265, 0.0]) < 1e-8)

  a.raw = '1 0 0\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.symmgroup == 'P 42/M N M'
  a.raw = a.raw
  assert a.symmgroup == 136

  a.raw = '1 1 0\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.ifhr == 1
  a.raw = a.raw
  assert a.ifhr == 1

  a.raw = '1 0 1\n'\
          'P 42/M N M\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert a.shift == 1
  a.raw = a.raw
  assert a.shift == 1

  a.raw = '1 0 2\n'\
          'P 42/M N M\n'\
          '5 5 5\n'\
          '4.63909875  2.97938395\n'\
          '2\n'\
          '22 0.0 0.0 0.0\n'\
          '8  3.061526467783E-01 3.061526467783E-01 0.0\n'
  assert all(abs(a.shift - 5/24.0) < 1e-8)
  a.raw = a.raw
  assert all(abs(a.shift - 5/24.0) < 1e-8)

  b = eval(repr(a), {'Crystal': Crystal, 'array': array})
  assert b.symmgroup == a.symmgroup
  assert all(abs(array(b.params) - a.params) < 1e-8)
  assert len(b) == len(a)
  assert len(b.atoms) == len(a.atoms)

def test_removeatoms():
  from numpy import all, arange
  from pylada.dftcrystal.geometry import RemoveAtoms
  a = RemoveAtoms()

  a.raw = '5\n1 2 3 4 5\n'
  assert all(a.labels == arange(5) + 1)
  a.raw = a.raw
  assert all(a.labels == arange(5) + 1)
  b = eval(repr(a), {'{0.__name__}'.format(RemoveAtoms): RemoveAtoms})
  assert all(b.labels == arange(5) + 1)

def test_modifysymmetry():
  from numpy import all, array
  from pylada.dftcrystal.geometry import ModifySymmetry
  a = ModifySymmetry()

  raw = '5\n1 A 2 A 3 0 4 0 5 0\n'
  a.raw = raw
  assert len(a.groups) == 2
  assert all(array(a.raw.split(), dtype='int32') == [5, 1, 1, 2, 1, 3, 2, 4, 2, 5, 2])
  a.raw = a.raw
  assert len(a.groups) == 2
  assert all(array(a.raw.split(), dtype='int32') == [5, 1, 1, 2, 1, 3, 2, 4, 2, 5, 2])
  b = eval(repr(a), {'{0.__name__}'.format(ModifySymmetry): ModifySymmetry})
  assert b.raw == a.raw


def test_slabinfo():
  from numpy import all, array
  from pylada.dftcrystal.geometry import Slabinfo

  a = Slabinfo()
  a.raw = '1 1 0'
  assert all(a.hkl == [1, 1, 0])
  assert all(array(a.raw.split(), dtype='int64') == [1, 1, 0])
  b = eval(repr(a), {'Slabinfo': Slabinfo})
  assert all(b.hkl == [1, 1, 0])

def test_slabcut():
  from numpy import all, array
  from pylada.dftcrystal.geometry import Slabcut

  a = Slabcut()
  a.raw = '1 1 0\n1 5\n'
  assert all(a.hkl == [1, 1, 0])
  assert a.isup == 1 and a.nl == 5
  assert all(array(a.raw.split(), dtype='int64') == [1, 1, 0, 1, 5])
  b = eval(repr(a), {'Slabcut': Slabcut})
  assert all(array(b.raw.split(), dtype='int64') == [1, 1, 0, 1, 5])

def test_displaceatoms():
  from numpy import all, abs
  from pylada.dftcrystal.geometry import DisplaceAtoms

  a = DisplaceAtoms()
  a.raw = '2\n1 0.1 0.1 0.1\n2 0.2 0.2 0.2'
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 1
  assert a.atoms[1].type == 2
  assert all(abs(a.atoms[0].pos - 0.1) < 1e-8)
  assert all(abs(a.atoms[1].pos - 0.2) < 1e-8)
  a.raw = a.raw
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 1
  assert a.atoms[1].type == 2
  assert all(abs(a.atoms[0].pos - 0.1) < 1e-8)
  assert all(abs(a.atoms[1].pos - 0.2) < 1e-8)

def test_insertatoms():
  from numpy import all, abs
  from pylada.dftcrystal.geometry import InsertAtoms

  a = InsertAtoms()
  a.raw = '2\n1 0.1 0.1 0.1\n2 0.2 0.2 0.2'
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'H'
  assert a.atoms[1].type == 'He'
  assert all(abs(a.atoms[0].pos - 0.1) < 1e-8)
  assert all(abs(a.atoms[1].pos - 0.2) < 1e-8)
  a.raw = a.raw
  assert len(a.atoms) == 2
  assert a.atoms[0].type == 'H'
  assert a.atoms[1].type == 'He'
  assert all(abs(a.atoms[0].pos - 0.1) < 1e-8)
  assert all(abs(a.atoms[1].pos - 0.2) < 1e-8)

def test_affinetransform():
  from numpy import all, abs, array
  from quantities import degree, angstrom
  from pylada.dftcrystal.geometry import AffineTransform
  a = AffineTransform()

  raw = '4\n'\
        '19 20 21 22\n'\
        '999 1\n'\
        '23 24\n'\
        '308\n'
  a.raw = raw
  assert a.vectrans is None
  assert a.bondtrans is None
  assert a.origin is None
  assert a.euler is None
  assert a.bondtoz is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert a.bondrot[0] == 23 and a.bondrot[1] == 24 and abs(a.bondrot[2].magnitude - 308) < 1e-8
  assert all(array(a.raw.split(), dtype='int32') == array(raw.split(), dtype='int32'))
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '999 1\n'\
        '23 24\n'\
        '0\n'
  a.raw = raw
  assert a.vectrans is None
  assert a.bondtrans is None
  assert a.origin is None
  assert a.euler is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert a.bondtoz[0] == 23 and a.bondtoz[1] == 24 
  assert all(array(a.raw.split(), dtype='int32') == array(raw.split(), dtype='int32'))
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '999 -1\n'\
        '23 24 25 1\n'
  a.raw = raw
  assert a.vectrans is None
  assert a.bondtrans is None
  assert a.origin is None
  assert a.bondtoz is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert a.euler[0] == 23 * degree and a.euler[1] == 24 * degree \
         and a.euler[2] == 25 * degree and a.euler[3] == 1
  assert a.raw == raw
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '1 999\n'
  a.raw = raw
  assert a.vectrans is None
  assert a.bondtrans is None
  assert a.euler is None
  assert a.bondtoz is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert a.origin == 1
  assert a.raw == raw
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '0 999\n'\
        '1 5 0.15\n'
  a.raw = raw
  assert a.vectrans is None
  assert a.origin is None
  assert a.euler is None
  assert a.bondtoz is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert a.bondtrans[0] == 1 and a.bondtrans[1] == 5 \
         and abs(a.bondtrans[2] - 0.15*angstrom) < 1e-8
  assert a.raw == raw
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '-1 999\n'\
        '1.0 5.0 0.15\n'
  a.raw = raw
  assert a.bondtrans is None
  assert a.origin is None
  assert a.euler is None
  assert a.bondtoz is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert all(abs(a.vectrans - [1, 5, 0.15] * angstrom) < 1e-8)
  assert a.raw == raw
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

  raw = '4\n'\
        '19 20 21 22\n'\
        '-1 -1\n'\
        '1.0 5.0 0.15\n'\
        '23 24 25 1\n'
  a.raw = raw
  assert a.bondtrans is None
  assert a.origin is None
  assert a.bondtoz is None
  assert a.bondrot is None
  assert all(array(a.labels) == [19, 20, 21, 22])
  assert all(abs(a.vectrans - [1, 5, 0.15] * angstrom) < 1e-8)
  assert a.euler[0] == 23 * degree and a.euler[1] == 24 * degree \
         and a.euler[2] == 25 * degree and a.euler[3] == 1
  assert a.raw == raw
  b = eval(repr(a), { 'AffineTransform': AffineTransform, 
                      'deg': degree, 'angstrom': angstrom, 'array': array})
  assert a.raw == b.raw

if __name__ == "__main__":
  test_crystal()
  test_removeatoms()
  test_modifysymmetry()
  test_slabinfo()
  test_slabcut()
  test_displaceatoms()
  test_insertatoms()
  test_affinetransform()

