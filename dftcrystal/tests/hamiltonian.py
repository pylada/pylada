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

def test():
  from collections import namedtuple
  from pylada.dftcrystal.hamiltonian import Dft
  from pylada.dftcrystal.input import print_input
  from pylada.dftcrystal.parse import parse
  Crystal = namedtuple('Crystal', ['dft'])
  
  a = Dft()
  crystal = Crystal(a)
  assert a.exchange is None
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is None
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert a.output_map(crystal=crystal) is None

  try: a.exchange = 'none'
  except: pass
  else: raise Exception()
  try: a.correlat = 'none'
  except: pass
  else: raise Exception()
  try: a.hybrid = 'none'
  except: pass
  else: raise Exception()
  a.spin = 'none'
  assert a.spin == False

  a.exchange = 'lda'
  crystal = Crystal(a)
  assert a.exchange == 'lda'
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is False
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert len(a.output_map(crystal=crystal)) == 1
  assert 'DFT' in a.output_map(crystal=crystal)
  assert len(a.output_map(crystal=crystal)['DFT']) == 1
  assert 'EXCHANGE' in a.output_map(crystal=crystal)['DFT']
  assert a.output_map(crystal=crystal)['DFT']['EXCHANGE'] == 'LDA'

  b = parse("DFT\nEXCHANGE\nPBE\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'pbe'
  assert a.correlat is None
  assert a.hybrid is None
  assert a.nonlocal.correlation is None
  assert a.nonlocal.exchange is None
  assert a.spin is None
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc

  a.b3lyp = True
  crystal = Crystal(a)
  assert a.exchange == 'becke'
  assert a.correlat == 'lyp'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is None
  assert a.b3lyp
  assert not a.b3pw
  assert not a.pbe0
  assert not a.soggaxc
  assert len(a.output_map(crystal=crystal)) == 1
  assert 'DFT' in a.output_map(crystal=crystal)
  assert len(a.output_map(crystal=crystal)['DFT']) == 1
  assert 'B3LYP' in a.output_map(crystal=crystal)['DFT']
  assert a.output_map(crystal=crystal)['DFT']['B3LYP'] is None

  b = parse("DFT\nB3LYP\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.b3lyp
  b = parse("DFT\nPBE0\nEND DFT\n", needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  crystal = Crystal(a)
  assert a.pbe0
  assert a.exchange == 'pbe'
  assert a.correlat == 'pbe'
  assert a.hybrid is None
  assert a.nonlocal.exchange is None
  assert a.nonlocal.correlation is None
  assert a.spin is None
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert len(a.output_map(crystal=crystal)) == 1
  assert 'DFT' in a.output_map(crystal=crystal)
  assert len(a.output_map(crystal=crystal)['DFT']) == 1
  assert 'PBE0' in a.output_map(crystal=crystal)['DFT']
  assert a.output_map(crystal=crystal)['DFT']['PBE0'] is None

  a.b3pw = True
  a.exchange = 'lda'
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is None
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  string = "DFT\nCORRELAT\nPWGGA\nEXCHANGE\nLDA\n"\
           "HYBRID\n" + str(float(20)) + "\n"     \
           "NONLOCAL\n" + str(float(0.9)) + " "   \
           + str(float(0.81)) + "\n" + "END DFT\n"
  assert print_input(a.output_map(crystal=Crystal(a))) == string
  assert print_input(a.output_map(crystal=Crystal(a))) == string

  b = parse(string, needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is None
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert print_input(a.output_map(crystal=Crystal(a))) == string

  string = "DFT\nCORRELAT\nPWGGA\nEXCHANGE\nLDA\n"\
           "HYBRID\n" + str(float(20)) + "\n"     \
           "NONLOCAL\n" + str(float(0.9)) + " "   \
           + str(float(0.81)) + "\nHELLO\nWORLD\n1\n" + "END DFT\n"
  b = parse(string, needstarters=False)
  a = Dft()
  a.read_input(b['DFT'])
  assert a.exchange == 'lda'
  assert a.correlat == 'pwgga'
  assert abs(a.hybrid - 20) < 1e-8
  assert abs(a.nonlocal.exchange - 0.9) < 1e-8
  assert abs(a.nonlocal.correlation - 0.81) < 1e-8
  assert a.spin is None
  assert not a.pbe0
  assert not a.b3lyp
  assert not a.b3pw
  assert not a.soggaxc
  assert 'hello' in a._input
  assert type(a._input['hello']) is bool and a._input['hello'] == True
  assert 'world' in a._input
  assert type(a._input['world']) is int and a._input['world'] == 1
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert 'HELLO' in string
  assert 'WORLD' in string 
  assert string[string.index('WORLD')+1] == '1'
  a.hello = False
  a.world = 2
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert 'HELLO' not in string
  assert 'WORLD' in string 
  assert string[string.index('WORLD')+1] == '2'

def test_grids():
  from numpy import array, all, abs
  from collections import namedtuple
  from pylada.dftcrystal.hamiltonian import Dft
  from pylada.dftcrystal.input import print_input
  Crystal = namedtuple('Crystal', ['dft'])
  a = Dft()
  assert a.radial.intervals is None
  assert a.radial.nbpoints is None
  assert a.angular.intervals is None
  assert a.angular.levels is None
  assert a.output_map(crystal=Crystal(a)) is None
  assert not a.lgrid
  assert not a.xlgrid
  assert not a.xxlgrid

  a.lgrid = True
  assert all(abs(array(a.angular.intervals) - [0.1667, 0.5, 0.9, 3.05, 9999.0]) < 1e-8)
  assert all(abs(array(a.angular.levels) - [2, 6, 8, 13, 8]) < 1e-8)
  assert all(abs(array(a.radial.intervals) - [4]) < 1e-8)
  assert all(abs(array(a.radial.nbpoints) - [75]) < 1e-8)
  assert print_input(a.output_map(crystal=Crystal(a))) == 'DFT\nLGRID\nEND DFT\n'


if __name__ == '__main__':
  test()
  test_grids()
