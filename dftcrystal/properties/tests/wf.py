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
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import join, exists
  from os import mkdir
  from pylada.dftcrystal import Crystal, Functional, Shell
  from pylada.dftcrystal.properties import Properties
  from pylada import default_comm

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
  functional.dft.spin = True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  directory =  mkdtemp() 
  if directory == '/tmp/test/':
    rmtree(directory)
    mkdir(directory)
  firstdir, seconddir = join(directory, '0'), join(directory, '1')
  try: 
     si = functional(crystal, outdir=firstdir, comm=default_comm)
     assert exists(join(firstdir, 'crystal.f9'))
     assert not exists(join(directory, 'crystal.f98'))
     properties = Properties(si)
     properties.rdfmwf = True
     result = properties(outdir=directory, workdir=join(directory, 'try0'))
     assert result.success
     assert exists(join(directory, 'crystal.f9'))
     assert not exists(join(directory, 'crystal.f98'))
     properties = Properties(directory)
     properties.fmwf = True
     result = properties(outdir=directory, workdir=join(directory, 'try1'))
     assert not exists(join(directory, 'crystal.f98'))
     result = properties(outdir=directory, workdir=join(directory, 'try1'), overwrite=True)
     assert result.success
     assert exists(join(directory, 'crystal.f98'))

  finally: 
    if directory != '/tmp/test/': rmtree(directory)
if __name__ == '__main__':
  test()
