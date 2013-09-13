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
  from pylada.dftcrystal.optgeom import OptGeom
  from pylada.dftcrystal.input import print_input
  Crystal = namedtuple('Crystal', ['optgeom'])

  a = OptGeom()
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = None
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is False
  assert a.itatocell is False
  assert a.intredun is False
  a.cellonly = True
  assert a.maxcycle is None
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.intredun is False
  a.maxcycle = 10
  assert a.maxcycle == 10
  assert a.fulloptg is False
  assert a.cellonly is True
  assert a.itatocell is False
  assert a.intredun is False

  assert a.enabled is False
  assert a.output_map(crystal=Crystal(a)) is None
  a.enabled = True
  assert a.enabled is True
  string = print_input(a.output_map(crystal=Crystal(a)))
  assert len(string.split('\n')) == 6
  string = string.split('\n')
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string

  a.cvolopt = True
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 7
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'CELLONLY' in string
  assert 'CVOLOPT' in string

  a.itatocell = True
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 6
  assert string[0] == 'OPTGEOM'
  assert string[-1] == ''
  assert string[-2] == 'END OPTGEOM'
  string = print_input(a.output_map(crystal=Crystal(a))).split()
  assert 'MAXCYCLE' in string
  assert string[string.index('MAXCYCLE')+1] == str(10)
  assert 'ITATOCELL' in string

  a.itatocell = False
  assert a.output_map(crystal=Crystal(a)) is not None
  string = print_input(a.output_map(crystal=Crystal(a))).split('\n')
  assert len(string) == 5
  assert string[0] == 'OPTGEOM'
  assert string[1] == 'MAXCYCLE'
  assert string[2] == '10'
  assert string[-2] == 'END OPTGEOM'
  assert string[-1] == ''

  a.enabled = False
  assert a.output_map(crystal=Crystal(a)) is None
 

def test_breaksym():
  from pylada.dftcrystal import Crystal, Functional, Shell, DisplaceAtoms, read

  functional = Functional()
  functional.basis['Si'] = [
      Shell('s', a0=[16120.0, 0.001959],
                 a1=[2426.0, 0.01493], 
                 a2=[553.9, 0.07285],
                 a3=[156.3, 0.2461], 
                 a4=[50.07, 0.4859],
                 a5=[17.02, 0.325]),
      Shell('sp', a0=[292.7, -0.002781, 0.004438],
                  a1=[69.87, -0.03571, 0.03267],
                  a2=[22.34, -0.115, 0.1347],
                  a3=[8.15, 0.09356, 0.3287],
                  a4=[3.135, 0.603, 0.4496]), 
      Shell('sp', 4.0, a0=[1.22, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.55, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.27, 1.0, 1.0]) ]

  functional.dft.pbe0 = True
  functional.fmixing = 30
  functional.shrink = 8, 16
  functional.levshift = 5, True
  functional.maxcycle = 600
  functional.optgeom.keepsymm = True
  functional.optgeom.enabled = True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  crystal.append('breaksym')
  crystal.append(DisplaceAtoms().add_atom(0.01, 0, 0, 1))
  string = functional.print_input(structure=crystal)
  assert string.splitlines()[10] == 'KEEPSYMM'
  cry, func = read(string)
  assert cry.is_breaksym
  assert func.optgeom.enabled 
  assert func.optgeom.keepsymm 

  functional.optgeom.keepsymm = False
  string = functional.print_input(structure=crystal)
  assert string.splitlines()[10] == 'OPTGEOM'
  cry, func = read(string)
  assert cry.is_breaksym
  assert func.optgeom.enabled 
  assert func.optgeom.breaksym 

  crystal[0] = 'keepsymm'
  string = functional.print_input(structure=crystal)
  assert string.splitlines()[11] == 'BREAKSYM'
  cry, func = read(string)
  assert not cry.is_breaksym
  assert func.optgeom.enabled 
  assert func.optgeom.breaksym 

  functional.optgeom.keepsymm = True
  string = functional.print_input(structure=crystal)
  assert string.splitlines()[11] == 'OPTGEOM'
  cry, func = read(string)
  assert not cry.is_breaksym
  assert func.optgeom.enabled 
  assert func.optgeom.keepsymm 

  string = crystal.print_input()
  crystal.append('keepsymm')
  crystal.append('keepsymm')
  assert crystal.print_input() == string
  crystal.append('breaksym')
  assert crystal.print_input() != string
  assert len(crystal.print_input().splitlines()) == len(string.splitlines()) + 1

if __name__ == "__main__":
  test()
  test_breaksym()

