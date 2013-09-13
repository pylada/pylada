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

string = """rutile
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
BREAKSYM
ATOMDISP
18
 1     0.0000000000     0.0000000000   -0.0720009980
 2     0.0000000000     0.0000000000   -0.1340625190
 3     0.0000000000     0.0000000000    0.1265375690
 4     0.0000000000    -0.0638880560    0.1707429180
 5     0.0000000000     0.0638880560    0.1707429180
 6     0.0000000000     0.0000000000   -0.0645623720
 7     0.0000000000     0.0000000000    0.0196680890
 8     0.0000000000     0.0000000000    0.0000000000
 9     0.0000000000     0.0000000000    0.0000000000
10     0.0000000000     0.0706330520    0.0000000000
11     0.0000000000    -0.0706330520    0.0000000000
12     0.0000000000     0.0000000000   -0.0196680890
13     0.0000000000     0.0000000000    0.0645623720
14     0.0000000000     0.0000000000    0.1340625190
15     0.0000000000     0.0000000000   -0.1265375690
16     0.0000000000    -0.0638880560   -0.1707429180
17     0.0000000000     0.0638880560   -0.1707429180
18     0.0000000000     0.0000000000    0.0720009980
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
99 0
END
DFT
B3LYP
XXLGRID
END
TOLINTEG
7 7 7 7 14
SHRINK
8 8
LEVSHIFT
5 1
FMIXING
60
PPAN
TOLDEE
7
EXCHSIZE
6937578
BIPOSIZE
9701800
SCFDIR
MAXCYCLE
300
END
"""
def test():
  """ Tests functional prints and reads itself """
  from pylada.dftcrystal.functional import Functional
  from pylada.dftcrystal import Crystal
  from pylada.dftcrystal.parse import parse
  parsed = parse(string)
  structure = Crystal()
  structure.read_input(parsed['rutile']['CRYSTAL'])
  a = Functional()
  a.read_input(parsed)
  assert a.scfdir 
  assert a.maxcycle == 300
  assert a.exchsize == 6937578
  # need structure otherwise parse can't find beginning of input.
  otherstring = a.print_input(structure=structure)
  otherparsed = parse(otherstring)
  b = Functional()
  b.read_input(otherparsed)
  assert otherstring == b.print_input(structure=structure)
  

def test_setprint():
  from pickle import loads, dumps
  from numpy import array, all
  from pylada.dftcrystal.input import SetPrint
  a = SetPrint()
  assert a.output_map() is None
  assert len(eval(repr(a), {'SetPrint': SetPrint}).options) == 0
  assert len(loads(dumps(a)).options) == 0
  b = SetPrint()
  b.read_input(a.output_map())
  assert len(b.options) == 0

  a[0] = 5
  assert len(a.options) == 1
  assert a.options[0] == 5
  assert a[0] == 5
  assert all(array(a.output_map()['setprint'].split(), dtype='int64') == [1, 0, 5])
  a[5] = 0
  assert len(a.options) == 2
  assert a.options[5] == 0
  assert a[5] == 0
  assert all(array(a.output_map()['setprint'].split(), dtype='int64') == [2, 0, 5, 5, 0])

  assert len(eval(repr(a), {'SetPrint': SetPrint}).options) == 2
  assert eval(repr(a), {'SetPrint': SetPrint}).options[0] == 5
  assert eval(repr(a), {'SetPrint': SetPrint}).options[5] == 0
  assert len(loads(dumps(a)).options) == 2
  assert loads(dumps(a)).options[0] == 5
  assert loads(dumps(a)).options[5] == 0
  b = SetPrint()
  b.read_input(a.output_map()['setprint'])
  assert len(b.options) == 2
  assert b[0] == 5 and b[5] == 0


if __name__ == '__main__': 
  test()
  test_setprint()
