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

def test_shrink():
  from pickle import loads, dumps
  from numpy import array, all
  from pylada.dftcrystal.electronic import Shrink
  from pylada.dftcrystal.molecule import Molecule

  a = Shrink()
  assert a.mp is None
  assert a.gallat is None
  assert a.output_map() is None
  assert eval(repr(a), {'Shrink': Shrink}).mp  is None
  assert eval(repr(a), {'Shrink': Shrink}).gallat  is None
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 1
  assert a.mp == 1 
  assert a.gallat is None
  assert 'shrink' in a.output_map()
  assert a.output_map()['shrink'] == a.raw
  assert a.output_map()['shrink'] == '1 1\n'
  assert eval(repr(a), {'Shrink': Shrink}).mp == 1
  assert eval(repr(a), {'Shrink': Shrink}).gallat  is None
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5
  assert a.mp == 5
  assert a.gallat is None
  assert a.output_map()['shrink'] == '5 5\n'
  assert repr(a) == 'Shrink(5)'
  assert repr(loads(dumps(a))) == repr(a)

  a.gallat = 10
  assert a.mp == 5
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '5 10\n'
  assert repr(a) == 'Shrink(5, 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 
  assert all(array(a.mp) == [5])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 1 1\n'
  assert repr(a) == 'Shrink([5], 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 6
  assert all(array(a.mp) == [5, 6])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 6 1\n'
  assert repr(a) == 'Shrink([5, 6], 10)'
  assert repr(loads(dumps(a))) == repr(a)

  a.mp = 5, 6, 7
  assert all(array(a.mp) == [5, 6, 7])
  assert a.gallat == 10
  assert a.output_map()['shrink'] == '0 10\n5 6 7\n'
  assert repr(a) == 'Shrink([5, 6, 7], 10)'
  assert repr(loads(dumps(a))) == repr(a)


  a.raw = '5 5\n'
  assert a.mp == 5 and a.gallat is None
  a.raw = '5 10\n'
  assert a.mp == 5 and a.gallat == 10
  a.raw = '0 5\n5 1 1'
  assert all(array(a.mp) == [5]) and a.gallat is None
  a.raw = '0 10\n5 2 1'
  assert all(array(a.mp) == [5, 2]) and a.gallat == 10
  a.raw = '0 10\n5 2 10'
  assert all(array(a.mp) == [5, 2, 10]) and a.gallat == 10
  a.raw = '0 10\n1 1 1'
  assert all(array(a.mp) == 1) and a.gallat == 10
  assert eval(repr(a), {'Shrink': Shrink}).mp  == 1
  assert eval(repr(a), {'Shrink': Shrink}).gallat == 10

  try: a.gallat = 'a'
  except ValueError: pass
  else: raise Exception()
  try: a.mp = 'a'
  except ValueError: pass
  else: raise Exception()
  try: a.mp = range(4)
  except ValueError: pass
  else: raise Exception()

  assert a.output_map(structure=Molecule()) is None

def test_levshift():
  from pickle import loads, dumps
  from quantities import hartree, UnitQuantity, eV, kbar
  from pylada.dftcrystal.electronic import LevShift, Electronic
  from pylada.error import ValueError

  a = LevShift()
  assert a.output_map() is None
  assert a.units == UnitQuantity('decihartree', 0.1*hartree)
  assert a.lock is None and a.shift is None
  assert repr(a) == 'LevShift()'
  assert repr(loads(dumps(a))) == repr(a)

  a.shift = 0.005 * eV
  d = {'LevShift': LevShift}
  assert a.output_map() is None
  assert a.lock is None
  assert abs(a.shift - 0.005 * eV) < 1e-8
  assert abs(a.shift.magnitude - 0.005 * 1./0.1 * eV.rescale(hartree).magnitude) < 1e-8
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)
  
  a.shift = 2
  assert a.output_map() is None
  assert a.lock is None
  assert abs(a.shift - 0.2 * hartree) < 1e-8
  assert abs(a.shift.magnitude - 2) < 1e-8
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a.lock = True
  assert a.raw == str(2) + ' ' + str(1)
  assert a.lock == True
  assert 'levshift' in a.output_map()
  assert a.output_map()['levshift'] == a.raw
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a.lock = False
  assert a.raw == str(2) + ' ' + str(0)
  assert abs(eval(repr(a), d).shift - a.shift) < 1e-8
  assert eval(repr(a), d).lock == a.lock
  assert repr(loads(dumps(a))) == repr(a)

  a = Electronic()
  assert a.levshift.raw == ''
  a.levshift = 2, False
  assert a.levshift.raw == str(2) + ' ' + str(0)
  assert abs(a.levshift[0] - 0.2 * hartree) < 1e-8
  assert a.levshift[1] is False

  try: a.levshift = 1
  except ValueError: pass
  else: raise Exception()
  try: a.levshift = 1*kbar, False
  except: pass
  else: raise Exception()

def test_spinlock():
  from collections import namedtuple
  from pickle import loads, dumps
  from pylada.dftcrystal.electronic import SpinLock, Electronic

  Crystal = namedtuple('Crystal', ['dft'])
  Dft = namedtuple('Dft', ['spin'])
  on, off = Crystal(Dft(True)), Crystal(Dft(False))
  a = Electronic()
  b = a._input['spinlock']
  assert a.spinlock.nspin is None
  assert a.spinlock.ncycles is None
  assert b.output_map() is None
  assert eval(repr(b), {'SpinLock': SpinLock}).output_map() is None
  assert eval(repr(b), {'SpinLock': SpinLock}).nspin is None
  assert eval(repr(b), {'SpinLock': SpinLock}).ncycles is None
  assert repr(loads(dumps(b))) == repr(b)

  a.spinlock = 5, 30
  assert a.spinlock[0] == 5
  assert a.spinlock[1] == 30
  assert a.spinlock.nspin == 5
  assert a.spinlock.ncycles == 30
  assert b.output_map(crystal=off) is None
  assert len(b.output_map(crystal=on)) == 1
  assert b.output_map(crystal=on).get('spinlock', 'a a') == '5 30' 
  assert eval(repr(b), {'SpinLock': SpinLock}).output_map(crystal=on)['spinlock'] == '5 30'
  assert eval(repr(b), {'SpinLock': SpinLock}).output_map(crystal=off) is None
  assert eval(repr(b), {'SpinLock': SpinLock}).nspin == 5
  assert eval(repr(b), {'SpinLock': SpinLock}).ncycles == 30
  assert repr(loads(dumps(b))) == repr(b)

def test_atomspin():
  from pickle import loads, dumps
  from collections import namedtuple
  from pylada.dftcrystal.electronic import AtomSpin, Electronic
  from pylada.error import TypeError

  class GuessP:
    def __init__(self, yes=False):
      self.yes = yes
    def output_map(self, **kwargs):
      return {'guessp': True} if self.yes else None
  class Scf:
    def __init__(self, yes=False):
      self._input = {'guessp': GuessP(yes)}

  Crystal = namedtuple('Crystal', ['dft', 'scf'])
  Dft = namedtuple('Dft', ['spin'])
  on, off = Crystal(Dft(True), Scf(False)), Crystal(Dft(False), Scf(False))
  a = Electronic()
  b = a._input['atomspin']
  assert a.atomspin.up is None
  assert a.atomspin.down is None
  assert a.atomspin.other is None
  assert b.output_map(crystal=on) is None
  assert b.output_map(crystal=off) is None

  a.atomspin.up = [1, 3, 5]
  assert a.atomspin.up == [1, 3, 5]
  assert b.output_map(crystal=off) is None
  assert 'atomspin' in b.output_map(crystal=on)
  c = b.output_map(crystal=on)['atomspin'].split()
  assert len(c) == 3*2+1
  assert c[0] == '3'
  assert all(u == '1' for u in c[2::2])
  assert c[1::2] == ['1', '3', '5']
  assert eval(repr(b), {'AtomSpin':AtomSpin}).output_map(crystal=on)           \
         == b.output_map(crystal=on)
  assert loads(dumps(b)).output_map(crystal=on)                                \
         == b.output_map(crystal=on)
  c = AtomSpin()
  c.read_input(b.output_map(crystal=on)['atomspin'])
  assert repr(c) == repr(b)

  a.atomspin.down = [2, 4, 6]
  assert a.atomspin.down == [2, 4, 6]
  assert b.output_map(crystal=off) is None
  assert 'atomspin' in b.output_map(crystal=on)
  c = b.output_map(crystal=on)['atomspin'].split()
  assert len(c) == 6*2+1
  assert c[0] == '6'
  assert all(u == '1' for u in c[2::2][:3])
  assert c[1::2][:3] == ['1', '3', '5']
  assert all(u == '-1' for u in c[2::2][3:])
  assert eval(repr(b), {'AtomSpin':AtomSpin}).output_map(crystal=on) \
         == b.output_map(crystal=on)
  assert loads(dumps(b)).output_map(crystal=on) \
         == b.output_map(crystal=on)
  assert c[1::2][3:] == ['2', '4', '6']
  c = AtomSpin()
  c.read_input(b.output_map(crystal=on)['atomspin'])
  assert repr(c) == repr(b)

  a.atomspin.other = [7, 8, 9]
  assert a.atomspin.other == [7, 8, 9]
  assert b.output_map(crystal=off) is None
  assert 'atomspin' in b.output_map(crystal=on)
  c = b.output_map(crystal=on)['atomspin'].split()
  assert len(c) == 9*2+1
  assert c[0] == '9'
  assert all(u == '1' for u in c[2::2][:3])
  assert c[1::2][:3] == ['1', '3', '5']
  assert all(u == '-1' for u in c[2::2][3:6])
  assert c[1::2][3:6] == ['2', '4', '6']
  assert all(u == '0' for u in c[2::2][6:9])
  assert c[1::2][6:9] == ['7', '8', '9']
  assert eval(repr(b), {'AtomSpin':AtomSpin}).output_map(crystal=on) \
         == b.output_map(crystal=on)
  assert loads(dumps(b)).output_map(crystal=on) \
         == b.output_map(crystal=on)
  c = AtomSpin()
  c.read_input(b.output_map(crystal=on)['atomspin'])
  assert repr(c) == repr(b)

  a.atomspin.up = None
  c = b.output_map(crystal=on)['atomspin'].split()
  assert len(c) == 6*2+1
  assert c[0] == '6'
  assert all(u == '-1' for u in c[2::2][:3])
  assert c[1::2][:3] == ['2', '4', '6']
  assert all(u == '0' for u in c[2::2][3:])
  assert c[1::2][3:] == ['7', '8', '9']
  assert eval(repr(b), {'AtomSpin':AtomSpin}).output_map(crystal=on) \
         == b.output_map(crystal=on)
  assert loads(dumps(b)).output_map(crystal=on) \
         == b.output_map(crystal=on)
  c = AtomSpin()
  c.read_input(b.output_map(crystal=on)['atomspin'])
  assert repr(c) == repr(b)
  
  a.atomspin.down = []
  c = b.output_map(crystal=on)['atomspin'].split()
  assert len(c) == 3*2+1
  assert c[0] == '3'
  assert all(u == '0' for u in c[2::2])
  assert c[1::2] == ['7', '8', '9']
  assert eval(repr(b), {'AtomSpin':AtomSpin}).output_map(crystal=on) \
         == b.output_map(crystal=on)
  assert loads(dumps(b)).output_map(crystal=on) \
         == b.output_map(crystal=on)
  c = AtomSpin()
  c.read_input(b.output_map(crystal=on)['atomspin'])
  assert c.output_map(crystal=on) == b.output_map(crystal=on)

  try: a.atomspin.down = 'a'
  except TypeError: pass
  else: raise Exception()
  try: a.atomspin.down = [0, 'a']
  except TypeError: pass
  else: raise Exception()
  try: a.atomspin.down = 5
  except TypeError: pass
  else: raise Exception()

def test_broyden():
  from pickles import loads, dumps
  from pylada.dftcrystal import Functional
  from pylada.dftcrystal.electronic import Broyden
  from pylada.error import ValueError, TypeError

  a = Functional()
  b = a.scf._input['broyden']
  assert a.broyden.w0 is None
  assert a.broyden.imix is None
  assert a.broyden.istart is None
  assert b.output_map() is None
  assert eval(repr(b), {'Broyden': Broyden}).w0 is None
  assert eval(repr(b), {'Broyden': Broyden}).imix is None
  assert eval(repr(b), {'Broyden': Broyden}).istart is None
  
  a.broyden.w0 = 1e-4
  assert abs(a.broyden.w0 - 1e-4) < 1e-8
  assert a.broyden.imix is None
  assert a.broyden.istart is None
  assert b.output_map() is None
  assert abs(eval(repr(b), {'Broyden': Broyden}).w0 - 1e-4) < 1e-8
  assert eval(repr(b), {'Broyden': Broyden}).imix is None
  assert eval(repr(b), {'Broyden': Broyden}).istart is None

  a.broyden.w0 = None
  a.broyden.imix = 50
  a.broyden.istart = 3
  assert a.broyden.w0 is None
  assert a.broyden.imix == 50
  assert a.broyden.istart == 3
  assert eval(repr(b), {'Broyden': Broyden}).w0 is None
  assert eval(repr(b), {'Broyden': Broyden}).imix == 50
  assert eval(repr(b), {'Broyden': Broyden}).istart == 3
  assert b.output_map() is None

  a.broyden.w0 = 1e-4
  a.broyden.imix = 50
  a.broyden.istart = 3
  assert abs(a.broyden.w0 - 1e-4) < 1e-8
  assert a.broyden.imix == 50
  assert a.broyden.istart == 3
  assert abs(eval(repr(b), {'Broyden': Broyden}).w0 - 1e-4) < 1e-8
  assert eval(repr(b), {'Broyden': Broyden}).imix == 50
  assert eval(repr(b), {'Broyden': Broyden}).istart == 3
  assert b.output_map() is not None
  assert 'broyden' in b.output_map() 
  assert len(b.output_map()['broyden'].split()) == 3
  assert abs(float(b.output_map()['broyden'].split()[0])-1e-4) < 1e-8
  assert int(b.output_map()['broyden'].split()[1]) == 50
  assert int(b.output_map()['broyden'].split()[2]) == 3

  a.broyden[0] = '1e-2'
  assert abs(a.broyden.w0 - 1e-2) < 1e-8
  assert abs(a.broyden[0] - 1e-2) < 1e-8
  a.broyden[1] = 25
  assert a.broyden.imix == 25
  assert a.broyden[1] == 25
  a.broyden[2] = 5
  assert a.broyden.istart == 5
  assert a.broyden[2] == 5

  assert abs(loads(dumps(b))[0] - 1e-2) < 1e-8
  assert loads(dumps(b))[1] == 25
  assert loads(dumps(b))[2] == 5

  a.broyden = None
  assert a.broyden.w0 is None
  assert a.broyden.imix is None
  assert a.broyden.istart is None
  a.broyden = ' 1e-2', 25, 5
  assert abs(a.broyden[0] - 1e-2) < 1e-8
  assert a.broyden[1] == 25
  assert a.broyden[2] == 5

  try: a.broyden = 0
  except TypeError: pass
  else: raise Exception()
  try: a.broyden = 0, 1
  except TypeError: pass
  else: raise Exception()
  try: a.broyden[1] = 0
  except ValueError: pass
  else: raise Exception()
  try: a.broyden[2] = 1
  except ValueError: pass
  else: raise Exception()


if __name__ == '__main__': 
  test_shrink()
  test_levshift()
  test_spinlock()
  test_atomspin()
  test_broyden()
