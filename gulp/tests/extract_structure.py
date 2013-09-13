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

def test_noopti_anatase(path):
  """ Test structure extraction. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom
  from pylada.gulp import Extract

  extract = Extract(join(path, 'noopti.gout'))
  assert not extract.optimize
  assert not extract.conp
  assert extract.qeq
  cell = array([[-1.897114,  1.897114,  4.884209],
                [ 1.897114, -1.897114,  4.884209],
                [ 1.897114,  1.897114, -4.884209]]) * angstrom
  assert all(abs(extract.cell - cell) < 1e-4)
  assert extract.cell is extract.input_cell
  assert all(abs(extract.structure.cell - cell.magnitude) < 1e-4)
  assert len(extract.structure) == 2
  assert extract.structure[0].type == 'Ti'
  assert all(abs(extract.structure[0].pos) < 1e-4)
  assert extract.structure[1].type == 'O'
  assert all(abs(extract.structure[1].pos - [0.99999295, 0.99999295, -0.99999295]) < 1e-4)
  assert extract.structure is extract.input_structure
  assert getattr(extract.structure[0], 'label', 0) == 1
  assert getattr(extract.structure[1], 'label', 0) == 2
  assert abs(getattr(extract.structure[0], 'occupancy', 0) - 1e0) < 1e-8
  assert abs(getattr(extract.structure[1], 'occupancy', 0) - 1e0) < 1e-8


def test_opti_anatase(path):
  """ Test structure optimization extraction. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom, e, eV
  from pylada.gulp import Extract

  extract = Extract(join(path, 'opti.gout'))
  assert extract.optimize
  assert extract.conp
  assert extract.qeq
  cell = array([[-1.897114,  1.897114,  4.884209],
                [ 1.897114, -1.897114,  4.884209],
                [ 1.897114,  1.897114, -4.884209]]) * angstrom
  assert all(abs(extract.input_cell - cell) < 1e-4)
  assert extract.cell is not extract.input_cell
  assert all(abs(extract.input_structure.cell - cell.magnitude) < 1e-4)
  assert len(extract.input_structure) == 2
  assert extract.input_structure[0].type == 'Ti'
  assert all(abs(extract.input_structure[0].pos) < 1e-4)
  assert extract.input_structure[1].type == 'O'
  assert all(abs(extract.input_structure[1].pos - [0.99999295, 0.99999295, -0.99999295]) < 1e-4)
  assert extract.structure is not extract.input_structure

  cell = array([[-1.924847,  1.924847,  4.53005 ],
                [ 1.924847, -1.924847,  4.53005 ],
                [ 1.924847,  1.924847, -4.53005 ]]) * angstrom
  assert all(abs(extract.cell - cell) < 1e-4)
  assert len(extract.structure) == 2
  assert extract.structure[0].type == 'Ti'
  assert all(abs(extract.structure[0].pos) < 1e-4)
  assert extract.structure[1].type == 'O'
  assert all(abs(extract.structure[1].pos - [0.98270828,  0.98270828, -0.98270828]) < 1e-4)
  assert getattr(extract.structure[0], 'label', 0) == 1
  assert getattr(extract.structure[1], 'label', 0) == 2
  assert abs(getattr(extract.structure[0], 'charge', 0) - 1.131064*e) < 1e-5
  assert abs(getattr(extract.structure[1], 'charge', 0) + 0.565532*e) < 1e-5
  assert abs(extract.structure.energy + 16.38287598 * eV) < 1e-5

def test_surface(path):
  """ Test extraction of surface structure. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom
  from pylada.gulp import Extract

  extract = Extract(join(path, 'surface.gout'))
  assert extract.optimize
  assert extract.conv
  assert extract.dimensionality == 2
  assert all(abs(extract.input_cell[:-1,:-1] - [[2.95796800, 0], [0, 6.48525389]]*angstrom) < 1e-4)
  assert all(abs(extract.input_cell[2] - [0, 0, 500.0]*angstrom) < 1e-8)
  assert all(abs(extract.input_cell[:, 2] - [0, 0, 500.0]*angstrom) < 1e-8)
  assert extract.cell is extract.input_cell
  input = array([a.pos for a in extract.input_structure])
  
  test  = array([[ 0.      ,  4.8639405 ,  6.262 ],
                 [ 0.      ,  1.6213135 ,  4.7179],
                 [ 1.478984,  4.8639405 ,  5.2297],
                 [ 1.478984,  2.83853082,  5.0944],
                 [ 1.478984,  0.40409618,  5.0944],
                 [ 0.      ,  4.8639405 ,  3.7261],
                 [ 0.      ,  1.6213135 ,  2.9489],
                 [ 0.      ,  4.8639405 ,  1.944 ],
                 [ 1.478984,  1.6213135 ,  1.4996],
                 [ 1.478984,  6.12938217,  1.6937],
                 [ 1.478984,  3.59849883,  1.6937],
                 [ 0.      ,  1.6213135 ,  0.4034],
                 [ 0.      ,  4.8639405 , -0.4034],
                 [ 0.      ,  1.6213135 , -1.944 ],
                 [ 1.478984,  4.8639405 , -1.4996],
                 [ 1.478984,  2.88675517, -1.6937],
                 [ 1.478984,  0.35587183, -1.6937],
                 [ 0.      ,  4.8639405 , -2.9489],
                 [ 0.      ,  1.6213135 , -3.7261],
                 [ 0.      ,  4.8639405 , -4.7179],
                 [ 1.478984,  1.6213135 , -5.2297],
                 [ 1.478984,  6.08115782, -5.0944],
                 [ 1.478984,  3.64672318, -5.0944],
                 [ 0.      ,  1.6213135 , -6.262 ]])
  assert all(abs(test-input) < 1e-6)
  input = array([a.pos for a in extract.structure])
  test  = array([[ 0.        ,  4.8639405 ,  6.742706  ],
                 [ 0.        ,  1.6213135 ,  5.039533  ],
                 [ 1.478984  ,  4.8639405 ,  5.510992  ],
                 [ 1.478984  ,  2.87609341,  5.344245  ],
                 [ 1.478984  ,  0.36653359,  5.344245  ],
                 [ 0.        ,  4.8639405 ,  3.959012  ],
                 [ 0.        ,  1.6213135 ,  3.091444  ],
                 [ 0.        ,  4.8639405 ,  2.03943   ],
                 [ 1.478984  ,  1.6213135 ,  1.660877  ],
                 [ 1.478984  ,  6.11509516,  1.740796  ],
                 [ 1.478984  ,  3.61278584,  1.740796  ],
                 [ 0.        ,  1.6213135 ,  0.403372  ],
                 [ 0.        ,  4.8639405 , -1.294272  ],
                 [ 0.        ,  1.6213135 , -2.930327  ],
                 [ 1.478984  ,  4.8639405 , -2.551775  ],
                 [ 1.478984  ,  2.87246816, -2.631687  ],
                 [ 1.478984  ,  0.37015884, -2.631687  ],
                 [ 0.        ,  4.8639405 , -3.982371  ],
                 [ 0.        ,  1.6213135 , -4.849909  ],
                 [ 0.        ,  4.8639405 , -5.930467  ],
                 [ 1.478984  ,  1.6213135 , -6.401891  ],
                 [ 1.478984  ,  6.1187269 , -6.235161  ],
                 [ 1.478984  ,  3.6091541 , -6.235161  ],
                 [ 0.        ,  1.6213135 , -7.6336    ]])

  assert all(abs(test-input) < 1e-6)


def test_optinosym(path):
  """ Test extraction of surface structure. """
  from os.path import join
  from numpy import array, abs, all
  from quantities import angstrom
  from pylada.gulp import Extract

  extract = Extract(join(path, 'opti_nosym.gout'))
  assert all( abs( extract.cell 
                   - array([[ 4.585767,  0.      , -0.      ],
                            [ 0.      ,  4.585767,  0.      ],
                            [-0.      ,  0.      ,  2.957967]]) * angstrom )
              < 1e-8 )
  input = array([a.pos for a in extract.structure])
  test = array([[ 0.        ,  0.        ,  0.        ],
                [ 2.2928835 ,  2.2928835 ,  1.4789835 ],
                [ 1.38293893,  1.38293893,  0.        ],
                [ 3.20282807,  3.20282807,  0.        ],
                [ 3.67582243,  0.90994457,  1.4789835 ],
                [ 0.90994457,  3.67582243,  1.4789835 ]])
  assert all(abs(test - input) < 1e-6)
  assert all( abs( extract.input_cell 
                   - array([[ 4.585767,  0.      , -0.      ],
                            [ 0.      ,  4.585767,  0.      ],
                            [-0.      ,  0.      ,  2.957968]]) * angstrom )
              < 1e-8 )
  input = array([a.pos for a in extract.input_structure])
  test = array([[ 0.        ,  0.        ,  0.        ],
                [ 2.2928835 ,  2.2928835 ,  1.478984  ],
                [ 1.40394632,  1.40394632,  0.        ],
                [ 3.18182068,  3.18182068,  0.        ],
                [ 3.69682982,  0.88893718,  1.478984  ],
                [ 0.88893718,  3.69682982,  1.478984  ]])
  assert all(abs(test - input) < 1e-6)


if __name__ == '__main__':
  from sys import argv
  from os.path import dirname, join
  test_noopti_anatase(join(dirname(argv[0]), 'data'))
  test_opti_anatase(join(dirname(argv[0]), 'data'))
  test_surface(join(dirname(argv[0]), 'data'))
  test_optinosym(join(dirname(argv[0]), 'data'))
