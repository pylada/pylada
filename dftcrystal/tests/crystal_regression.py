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

def test_doubleatomsymm():
  """ Bugs when calling atomsymm twice. """
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import exists
  from pylada.dftcrystal import Crystal
  from pylada.misc import Changedir


  c = Crystal(136, 4.63909875, 2.97938395, 
              ifhr=0, 
              shift=0)                                  \
             .add_atom(0, 0, 0, 'Ti')                   \
             .add_atom(0.306153, 0.306153, 0, 'O')      \
             .append('ATOMSYMM')
  
  directory = '/tmp/test' #mkdtemp()
  if directory == '/tmp/test':
    if exists(directory): rmtree(directory)
    with Changedir(directory) as cwd: pass
  try: 
      c.append('ATOMSYMM')
      c.eval()
  finally:
    if directory != '/tmp/test': rmtree(directory)

if __name__ == '__main__':
  test_doubleatomsymm()
