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

def test_general():
  from numpy import all, abs, array
  from quantities import angstrom
  from pylada.dftcrystal.basis import Shell
  from pylada.error import TypeError
  bohr = 0.5291772083*angstrom # crystal conversion factor
  invbohr = 1e0/bohr/bohr
  invangs = 1e0/angstrom/angstrom

  assert Shell(0).type == 's'  and abs(Shell(0).charge -  2.0) < 1e-8\
         and Shell(0).type_as_int == 0
  assert Shell(1).type == 'sp' and abs(Shell(1).charge -  8.0) < 1e-8\
         and Shell(1).type_as_int == 1
  assert Shell(2).type == 'p'  and abs(Shell(2).charge -  6.0) < 1e-8\
         and Shell(2).type_as_int == 2
  assert Shell(3).type == 'd'  and abs(Shell(3).charge - 10.0) < 1e-8\
         and Shell(3).type_as_int == 3
  assert Shell(4).type == 'f'  and abs(Shell(4).charge - 18.0) < 1e-8\
         and Shell(4).type_as_int == 4

  assert len(Shell(0).append(0.5, 0.1)) == 1
  assert len(Shell(0).append(0.5, 0.1).append(0.6, 0.1)) == 2
  assert abs(Shell(0).append(0.5, 0.1).functions[0][0] - 0.5*invbohr) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][0] - 0.5*invangs) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][0].magnitude - 0.140014258892) < 1e-8
  assert abs(Shell(0).append(0.5*invangs, 0.1).functions[0][1] - 0.1) < 1e-8
  assert len(Shell(0).append(0.5*invangs, 0.1).functions[0]) == 2
  assert len(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0]) == 3
  assert abs(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0][1] - 0.1) < 1e-8
  assert abs(Shell(1).append(0.5*invangs, 0.1, 0.2).functions[0][2] - 0.2) < 1e-8
  try: Shell(0).append(0.5*invangs, 0.1, 0.2)
  except TypeError: pass
  else: raise Exception()
  try: Shell(1).append(0.5*invangs, 0.1)
  except TypeError: pass
  else: raise Exception()

  a = Shell()
  raw = '0 0 8 2. 1.\n'\
        '225338.0 0.000228\n'\
        '32315.0 0.001929\n'\
        '6883.61 0.011100\n'\
        '1802.14 0.05\n'\
        '543.063 0.17010\n'\
        '187.549 0.369\n'\
        '73.2133 0.4033\n'\
        '30.3718 0.1445\n'
  a.raw = raw
  assert a.type == 's'
  assert abs(a.charge - 2) < 1e-8
  assert len(a) == 8
  assert all(array([len(u) for u in a.functions]) == 2)
  assert all(abs( array([float(u[0]) for u in a.functions])
                  - [225338.0, 32315.0, 6883.61, 1802.14, 543.063, 187.549, 73.2133, 30.3718] ) < 1e-8)
  assert all(abs( array([float(u[1]) for u in a.functions])
                  - [2.28e-4, 1.929e-3, 1.11e-2, 5e-2, 1.701e-1, 3.69e-1, 4.033e-1, 1.445e-1 ] ) < 1e-8)
  a_ = array(a.raw.split(), dtype='float64')
  assert all(abs(array(raw.split(), dtype='float64') - a_) < 1e-8)
  b = eval(repr(a), {'Shell': Shell})
  assert all(abs(array(b.raw.split(), dtype='float64') - a_) < 1e-8)

  raw = '0 1 6 8. 1.\n'\
        '554.042 -0.0059 0.0085\n'\
        '132.525 -0.0683 0.0603\n'\
        '43.6801 -0.1245 0.2124\n'\
        '17.2243 0.2532 0.3902\n'\
        '7.2248 0.6261 0.4097\n'\
        '2.4117 0.282 0.2181\n'
  a.raw = raw
  assert a.type == 'sp'
  assert abs(a.charge - 8) < 1e-8
  assert len(a) == 6
  assert all(array([len(u) for u in a.functions]) == 3)
  assert all(abs( array([float(u[0]) for u in a.functions])
                  - [554.042, 132.525, 43.6801, 17.2243, 7.2248, 2.4117] ) < 1e-8)
  assert all(abs( array([float(u[1]) for u in a.functions])
                  - [-0.0059, -0.0683, -0.1245, 0.2532, 0.6261, 0.282] ) < 1e-8)
  assert all(abs( array([float(u[2]) for u in a.functions])
                  - [0.0085, 0.0603, 0.2124, 0.3902, 0.4097, 0.2181] ) < 1e-8)
  a_ = array(a.raw.split(), dtype='float64')
  assert all(abs(array(raw.split(), dtype='float64') - a_) < 1e-8)
  b = eval(repr(a), {'Shell': Shell})
  assert all(abs(array(b.raw.split(), dtype='float64') - a_) < 1e-8)

input = """
CRYSTAL
0 0 0
136
4.63909875  2.97938395
2
22 0.0 0.0 0.0
8  3.061526467783E-01 3.061526467783E-01 0.0
SLAB
1 1 0
2 9
OPTGEOM
MAXCYCLE
10
END
END
22 7
0 0 8 2. 1.
225338.0 0.000228
32315.0 0.001929
6883.61 0.011100
1802.14 0.05
543.063 0.17010
187.549 0.369
73.2133 0.4033
30.3718 0.1445
0 1 6 8. 1.
554.042 -0.0059 0.0085
132.525 -0.0683 0.0603
43.6801 -0.1245 0.2124
17.2243 0.2532 0.3902
7.2248 0.6261 0.4097
2.4117 0.282 0.2181
0 1 4 8. 1.
24.4975 0.0175 -0.0207
11.4772 -0.2277 -0.0653
4.4653 -0.7946 0.1919
1.8904 1.0107 1.3778
0 1 1 2. 1.
0.8126 1.0 1.0
0 1 1 0. 1.
0.3297 1.0 1.0
0 3 4 2. 1.
16.2685 0.0675
4.3719 0.2934
1.4640 0.5658
0.5485 0.5450
0 3 1 0. 1.
0.26 1.0
8 5
0 0 8 2. 1.
8020.0 0.00108
1338.0 0.00804
255.4 0.05324
69.22 0.1681
23.90 0.3581
9.264 0.3855
3.851 0.1468
1.212 0.0728
0 1 4 6. 1.
49.43 -0.00883 0.00958
10.47 -0.0915 0.0696
3.235 -0.0402 0.2065
1.217 0.379 0.347
0 1 1 0. 1.
0.4567 1.0 1.0
0 1 1 0. 1.
0.1843 1.0 1.0
0 3 1 0. 1.
 0.6 1.0
1 3
0 0 3 1. 1.
18.7311370 0.03349460
2.82539370 0.234726950
0.640121700 0.813757330
0 0 1 0. 1.
0.16127780 1.0
0 2 1 0. 1.
1.1 1.0
99 0
BREAKSYM
GHOSTS
2
1 2
END
DFT
B3LYP
END
MAXCYCLE
300
END """

def test_read_input():
  from numpy import array, all
  from pylada.dftcrystal import Functional
  from pylada.dftcrystal.parse import parse
  a = Functional()
  parsed = parse(input)['']['BASISSET']
  a.basis.read_input(parsed)
  assert len(a.basis.ghosts) == 2
  assert a.basis.ghosts[0] == 1
  assert a.basis.ghosts[1] == 2
  assert a.basis.ghosts.breaksym is True
  assert len(a.basis) == 3
  assert set(a.basis) == set(['H', 'O', 'Ti'])
  assert len(a.basis['H']) == 3
  assert abs(a.basis['H'][0].charge - 1) < 1e-8
  assert len(a.basis['H'][0].functions) == 3
  assert a.basis['H'][0].type == 's'
  assert all(abs(array(a.basis['H'][0].functions[0]) - [18.7311370, 0.03349460]) < 1e-8)
  assert all(abs(array(a.basis['H'][0].functions[1]) - [2.82539370, 0.234726950]) < 1e-8)
  assert all(abs(array(a.basis['H'][0].functions[2]) - [0.640121700, 0.813757330]) < 1e-8)

def test_output_map():
  from pylada.dftcrystal import Functional, Crystal
  from pylada.dftcrystal.parse import parse

  a = Functional()
  parsed = parse(input)['']
  a.basis.read_input(parsed['BASISSET'])
  structure = Crystal()
  structure.read_input(parsed['CRYSTAL'])
  o = a.basis.output_map(structure=structure)
  assert len(o) == 1
  assert o[0] == ('ghosts', '2\n1 2')
  a.basis.ghosts.breaksym = False
  o = a.basis.output_map(structure=structure)
  assert len(o) == 2
  assert o[0] == ('keepsymm', True)
  assert o[1] == ('ghosts', '2\n1 2')
  b = Functional()
  b.basis.raw = o.prefix
  assert set(b.basis) == set(['O', 'Ti'])
  assert repr(b.basis['O']) == repr(a.basis['O'])
  assert repr(b.basis['Ti']) == repr(a.basis['Ti'])

if __name__ == "__main__":
  test_general()
  test_read_input()
  test_output_map()
