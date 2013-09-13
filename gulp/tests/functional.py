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

def test_input():
  from numpy import array, all, abs
  from pylada.error import ValueError
  from pylada.gulp import Functional

  functional = Functional()
  functional.morse['Ti', 'O']  = [1.0279493, 3.640737, 1.88265, 0, 25]
  functional.morse['Ti', 'Ti'] = [0.00567139, 1.5543, 4.18784, 0, 25]
  functional.morse['O', 'O']   = [0.00567139, 1.5543, 4.18784, 0, 25] 

  # should fail since no interaction enabled yet.
  try: functional.input_string()
  except ValueError: pass
  else: raise Exception()

  functional.morse.enabled = True
  string = functional.input_string().splitlines()
  assert string[:2] == [''] * 2
  assert string[2] == 'morse'
  assert string[3].split()[:4] == ['O', 'core'] * 2
  assert all(abs(array( string[3].split()[4:], dtype='float64')
                        - functional.morse['O', 'O'] ) < 1e-3) 
  assert string[4].split()[:4] == ['Ti', 'core'] * 2
  assert all(abs(array( string[4].split()[4:], dtype='float64')
                        - functional.morse['Ti', 'Ti'] ) < 1e-3) 
  assert string[5].split()[:4] == ['O', 'core', 'Ti', 'core']
  assert all(abs(array( string[5].split()[4:], dtype='float64')
                        - functional.morse['Ti', 'O'] ) < 1e-3) 

def test_repr():
  """ Test representability. """
  from numpy import all, array, abs
  from pylada.gulp import Functional
  
  functional = Functional()
  b = {}; exec(repr(functional)) in b; b = b['functional']
  assert not b.morse.enabled
  assert not b.opti
  assert not b.conv
  assert not b.conp
  assert not b.isotropic
  assert not b.cellonly
  assert not b.buckingham.enabled
  assert len(b.morse) == 0
  assert len(b.buckingham) == 0

  functional.morse['Ti', 'O']  = [1.0279493, 3.640737, 1.88265, 0, 25]
  functional.morse['Ti', 'Ti'] = [0.00567139, 1.5543, 4.18784, 0, 25]
  functional.morse['O', 'O']   = [0.00567139, 1.5543, 4.18784, 0, 25] 
  b = {}; exec(repr(functional)) in b; b = b['functional']
  assert not b.morse.enabled
  assert not b.opti
  assert not b.conv
  assert not b.conp
  assert not b.isotropic
  assert not b.cellonly
  assert len(b.morse) == 3
  assert not b.buckingham.enabled
  assert len(b.buckingham) == 0
  assert all(abs(array( functional.morse['Ti', 'O'])  
                        - [1.0279493, 3.640737, 1.88265, 0, 25] ) < 1e-4)
  assert all(abs(array( functional.morse['Ti', 'Ti']) 
                        - [0.00567139, 1.5543, 4.18784, 0, 25] ) < 1e-4)
  assert all(abs(array( functional.morse['O', 'O'])
                        - [0.00567139, 1.5543, 4.18784, 0, 25] ) < 1e-4)

  functional.morse.enabled = True
  functional.opti          = True
  functional.cellonly      = True
  b = {}; exec(repr(functional)) in b; b = b['functional']
  assert b.morse.enabled
  assert b.opti
  assert not b.conv
  assert not b.conp
  assert not b.isotropic
  assert not b.buckingham.enabled
  assert len(b.buckingham) == 0
  assert b.cellonly
  assert len(b.morse) == 3
  assert all(abs(array( functional.morse['Ti', 'O'])  
                        - [1.0279493, 3.640737, 1.88265, 0, 25] ) < 1e-4)
  assert all(abs(array( functional.morse['Ti', 'Ti']) 
                        - [0.00567139, 1.5543, 4.18784, 0, 25] ) < 1e-4)
  assert all(abs(array( functional.morse['O', 'O'])
                        - [0.00567139, 1.5543, 4.18784, 0, 25] ) < 1e-4)

def test_filesetup():
  """ Tests file setup. """
  from tempfile import mkdtemp
  from os.path import exists, join, samefile, lexists
  from shutil import rmtree
  from pylada.gulp import Functional
  from pylada.crystal import Structure
  from pylada.tools import create_directory, get_section_from_file

  functional = Functional()
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
    iterator = functional.iter(structure, outdir=directory)
    assert exists(directory)
    assert not exists(join(directory, 'gulp.out'))
    assert not exists(join(directory, 'gulp.in'))
    assert not exists(join(directory, 'gulp.err'))
    
    program = iterator.next()
    assert exists(join(directory, 'gulp.out'))
    assert exists(join(directory, 'gulp.in'))
    assert not exists(join(directory, 'gulp.err'))
    with open(join(directory, 'gulp.out'), 'r') as file:
      get_section_from_file(file, 'input file')
      get_section_from_file(file, 'functional')

    rmtree(directory)
    create_directory(directory)

    workdir = join(directory, 'wtf')
    iterator = functional.iter( structure, outdir=directory,
                                workdir=workdir )
    assert exists(directory)
    assert not exists(join(directory, 'gulp.out'))
    assert not exists(join(directory, 'gulp.in'))
    assert not exists(join(directory, 'gulp.err'))

    program = iterator.next()
    assert exists(join(directory, 'gulp.out'))
    assert exists(join(directory, 'gulp.in'))
    assert not exists(join(directory, 'gulp.err'))
    with open(join(directory, 'gulp.out'), 'r') as file:
      get_section_from_file(file, 'input file')
      get_section_from_file(file, 'functional')
    assert samefile(workdir, join(directory, 'workdir'))
    assert exists(join(workdir, 'gulp.out'))
    assert samefile(join(workdir, 'gulp.out'), join(directory, 'gulp.out'))
    assert exists(join(workdir, 'gulp.in'))
    assert not samefile(join(workdir, 'gulp.in'), join(directory, 'gulp.in'))
    assert lexists(join(workdir, 'gulp.err'))
    assert exists(join(directory, '.pylada_is_running'))

    program.onfinish()
    assert not exists(join(directory, '.pylada_is_running'))
    assert exists(workdir)
    assert exists(join(directory, 'workdir'))
    assert samefile(join(workdir, 'gulp.out'), join(directory, 'gulp.out'))

  finally:
    if directory != '/tmp/test/': rmtree(directory)

if __name__ == '__main__':
  test_input()
  test_repr()
  test_filesetup()
