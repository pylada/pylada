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

def test_chemod():
  from pickle import loads, dumps
  from numpy import all, array, abs
  from pylada.dftcrystal import Crystal, Functional, Shell
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
  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')

  kwargs = {'crystal': functional, 'structure': crystal}
  a = functional.basis._input['chemod']
  assert len(functional.basis.chemod) == 0
  assert a.output_map(**kwargs) is None
  assert functional.basis.chemod.breaksym == True
  assert loads(dumps(a)).breaksym is True
  assert len(loads(dumps(a))) == 0
 
  functional.basis.chemod.keepsym = True
  functional.basis.chemod[1] = [2.0, 8.0, 4.0, 0.0, 0.0]
  assert 'chemod' in a.output_map(**kwargs)
  assert all( abs( array(a.output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  assert a.modisymm(crystal) is None
  assert loads(dumps(a)).breaksym is False
  assert all( abs( array(loads(dumps(a)).output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  functional.basis.chemod[2] = [2.0, 8.0, 4.0, 0.0, 0.0]
  assert 'chemod' in a.output_map(**kwargs)
  assert all( abs( array(a.output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  assert loads(dumps(a)).breaksym is False
  assert all( abs( array(loads(dumps(a)).output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  assert a.modisymm(crystal) is None
  del functional.basis.chemod[1]
  assert 'chemod' in a.output_map(**kwargs)
  assert all( abs( array(a.output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  assert loads(dumps(a)).breaksym is False
  assert all( abs( array(loads(dumps(a)).output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [2.0,  1.0, 2.0, 8.0, 4, 0, 0, 2, 2, 8, 4, 0, 0] ) < 1e-8 )
  assert a.modisymm(crystal) is None

  functional.basis.chemod.breaksym = True
  assert 'chemod' in a.output_map(**kwargs)
  assert all( abs( array(a.output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [1.0,  2.0, 2.0, 8.0, 4, 0, 0] ) < 1e-8 )
  assert loads(dumps(a)).breaksym is True
  assert all( abs( array(loads(dumps(a)).output_map(**kwargs)['chemod'].split(), dtype='float64')
                   - [1.0,  2.0, 2.0, 8.0, 4, 0, 0] ) < 1e-8 )
  assert a.modisymm(crystal) is not None
  assert [2] in a.modisymm(crystal).groups 
  assert [1] in a.modisymm(crystal).groups 
  assert len(a.modisymm(crystal).groups) == 2

if __name__ == '__main__': 
  test_chemod()
