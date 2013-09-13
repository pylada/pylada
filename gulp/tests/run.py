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
  """ Tests file setup. """
  from tempfile import mkdtemp
  from os.path import exists
  from shutil import rmtree
  from quantities import eV
  from pylada.gulp import Functional
  from pylada.crystal import Structure
  from pylada.tools import create_directory

  functional = Functional()
  functional.qeq = True
  functional.morse.enabled     = True
  functional.morse['Ti', 'O']  = [1.0279493, 3.640737, 1.88265, 0, 25]
  functional.morse['Ti', 'Ti'] = [0.00567139, 1.5543, 4.18784, 0, 25]
  functional.morse['O', 'O']   = [0.00567139, 1.5543, 4.18784, 0, 25] 

  structure = Structure( 4.6391, 0, 0,
                         0, 4.6391, 0,
                         0, 0, 2.97938,
                         scale=1, symmgroup=136 )                              \
                .add_atom(0, 0, 0, 'Ti', asymmetric=True)                      \
                .add_atom(2.31955, 2.31955, 1.48969, 'Ti', asymmetric=False)   \
                .add_atom(1.42027, 1.42027, 0, 'O', asymmetric=True)           \
                .add_atom(-1.42027, -1.42027, 0, 'O', asymmetric=False)        \
                .add_atom(-0.899275, 0.899275, 1.48969, 'O', asymmetric=False) \
                .add_atom(0.899275, -0.899275, 1.48969, 'O', asymmetric=False)

  directory = mkdtemp()
  if directory == '/tmp/test/':
    if exists(directory): rmtree(directory)
    create_directory(directory)
  try:
    result = functional(structure, directory)
    assert result._finished
    assert result.success
    assert not result.optimize
    assert not result.conp
    assert abs(result.energy + 9.84502039*eV) < 1e-5
  finally:
    if directory != '/tmp/test/': rmtree(directory)

if __name__ == '__main__':
  test()
