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

def test_supercell(path):
  from os.path import join
  from numpy import all, abs, array, set_printoptions
  from pylada.dftcrystal import read, Supercell, Crystal, Slabcut, BreakSym,    \
                              KeepSymm, DisplaceAtoms, InsertAtoms
  set_printoptions(precision=12)
  crystal, functional = read(join(path, 'supercell.d12'))

  assert isinstance(crystal, Crystal)
  assert all(abs(array(crystal.params)  - [4.63909875, 2.97938395]) < 1e-6)
  assert crystal.symmgroup == 136
  assert len(crystal) == 8

  assert isinstance(crystal[0], Slabcut)
  assert crystal[0].isup == 2
  assert crystal[0].nl == 15
  assert all(crystal[0].hkl == [1, 1, 0])

  assert isinstance(crystal[1], BreakSym)

  assert isinstance(crystal[2], DisplaceAtoms)
  assert all(array([u.type for u in crystal[2].atoms]) == range(1, 31))
  disp = array([[ 0.         ,  0.         ,  0.00496008 ],
                [ 0.         ,  0.         , -0.152840221],
                [ 0.         ,  0.         ,  0.216529621],
                [ 0.         , -0.052728927,  0.185587552],
                [ 0.         ,  0.052728927,  0.185587552],
                [ 0.         ,  0.         , -0.004768893],
                [ 0.         ,  0.         ,  0.034440171],
                [ 0.         ,  0.         ,  0.136085018],
                [ 0.         ,  0.         , -0.069420272],
                [ 0.         ,  0.035782986,  0.022806847],
                [ 0.         , -0.035782986,  0.022806847],
                [ 0.         ,  0.         ,  0.018900892],
                [ 0.         ,  0.         ,  0.002450591],
                [ 0.         ,  0.         ,  0.         ],
                [ 0.         ,  0.         ,  0.         ],
                [ 0.         , -0.027503191,  0.         ],
                [ 0.         ,  0.027503191,  0.         ],
                [ 0.         ,  0.         , -0.002450591],
                [ 0.         ,  0.         , -0.018900892],
                [ 0.         ,  0.         , -0.136085018],
                [ 0.         ,  0.         ,  0.069420272],
                [ 0.         ,  0.035782986, -0.022806847],
                [ 0.         , -0.035782986, -0.022806847],
                [ 0.         ,  0.         , -0.034440171],
                [ 0.         ,  0.         ,  0.004768893],
                [ 0.         ,  0.         ,  0.152840221],
                [ 0.         ,  0.         , -0.216529621],
                [ 0.         , -0.052728927, -0.185587552],
                [ 0.         ,  0.052728927, -0.185587552],
                [ 0.         ,  0.         , -0.00496008 ]])
  assert all(abs( array([u.pos for u in crystal[2].atoms]) - disp) < 1e-8)

  assert isinstance(crystal[3], BreakSym)

  assert isinstance(crystal[4], DisplaceAtoms)
  assert all(array([u.type for u in crystal[4].atoms]) == range(1, 31))
  disp = array([[ 0.         ,  0.         , -0.005827511],
                [ 0.         ,  0.         , -0.00382378 ],
                [ 0.         ,  0.         , -0.000865488],
                [ 0.         ,  0.001041455, -0.002644481],
                [ 0.         , -0.001041455, -0.002644481],
                [ 0.         ,  0.         , -0.004243966],
                [ 0.         ,  0.         ,  0.001098565],
                [ 0.         ,  0.         ,  0.00118787 ],
                [ 0.         ,  0.         , -0.002160436],
                [ 0.         ,  0.000613051, -0.000776595],
                [ 0.         , -0.000613051, -0.000776595],
                [ 0.         ,  0.         ,  0.00162428 ],
                [ 0.         ,  0.         , -0.000458678],
                [ 0.         ,  0.         ,  0.         ],
                [ 0.         ,  0.         ,  0.         ],
                [ 0.         , -0.000733281,  0.         ],
                [ 0.         ,  0.000733281,  0.         ],
                [ 0.         ,  0.         ,  0.000458678],
                [ 0.         ,  0.         , -0.00162428 ],
                [ 0.         ,  0.         , -0.00118787 ],
                [ 0.         ,  0.         ,  0.002160436],
                [ 0.         ,  0.000613051,  0.000776595],
                [ 0.         , -0.000613051,  0.000776595],
                [ 0.         ,  0.         , -0.001098565],
                [ 0.         ,  0.         ,  0.004243966],
                [ 0.         ,  0.         ,  0.00382378 ],
                [ 0.         ,  0.         ,  0.000865488],
                [ 0.         ,  0.001041455,  0.002644481],
                [ 0.         , -0.001041455,  0.002644481],
                [ 0.         ,  0.         ,  0.005827511]])
  assert all(abs( array([u.pos for u in crystal[4].atoms]) - disp) < 1e-8)

  assert crystal[5].keyword.lower() == 'supercel'
  assert isinstance(crystal[5], Supercell)
  assert all(abs(crystal[5].matrix - [[2, 0], [0, 1]]) < 1e-8)

  assert isinstance(crystal[6], KeepSymm)

  assert isinstance(crystal[7], InsertAtoms)
  assert len(crystal[7].atoms) == 2
  assert crystal[7].atoms[0].type == 'H'
  assert crystal[7].atoms[1].type == 'H'
  assert all(abs(crystal[7].atoms[0].pos -  [0.0, 3.28033818472, 8.731578688]) < 1e-8)
  assert all(abs(crystal[7].atoms[1].pos -  [2.97938395, 3.28033818472, 8.731578688]) < 1e-8)


  assert all(functional.atomspin.up == [3, 4, 15, 16, 39, 40, 51, 52])
  assert functional.biposize             == 9701800
  assert functional.dft.b3lyp            == True
  assert functional.dft.spin             == True
  assert functional.dft.xxlgrid          == True
  assert functional.exchsize             == 6937578
  assert functional.fmixing              == 60
  assert functional.levshift.lock        == True
  assert functional.levshift.shift       == 5.0
  assert functional.maxcycle             == 300
  assert functional.mpp                  == True
  assert functional.nobipola             == True
  assert functional.nofmwf               == True
  assert functional.optgeom.breaksym     == False
  assert functional.optgeom.enabled      == True
  assert functional.optgeom.maxcycle     == 10
  assert functional.ppan                 == True
  assert functional.scfdir               == True
  assert all(functional.shrink.mp == [4, 8])
  assert functional.shrink.gallat == 8
  assert functional.spinlock.ncycles     == 30
  assert functional.spinlock.nspin       == 4
  assert functional.title                == 'rutile'
  assert functional.toldee               == 7
  assert all(functional.tolinteg         == [7, 7, 7, 7, 14])

if __name__ == '__main__':
  from sys import argv
  from os.path import dirname, join
  test_supercell(join(dirname(argv[0]), 'data'))
