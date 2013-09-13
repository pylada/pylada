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
  from numpy import array, all, abs
  from shutil import rmtree
  from os.path import exists
  from os import mkdir
  from pylada.dftcrystal import Crystal, relax, Shell, DisplaceAtoms
  from pylada import default_comm

  functional = relax.Relax()
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
  functional.optgeom.fulloptg = True
  functional.shrink = 8, 16
  functional.levshift = 5, True
  functional.maxcycle = 600
  functional.guessp = True
  functional.optgeom.keepsymm = True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  crystal.append('fraction')
  crystal.append('keepsymm')
  crystal.append(DisplaceAtoms().add_atom(0.01, 0, 0, 1))
  directory = mkdtemp()
  if directory == '/tmp/test/' and exists(directory):
    rmtree(directory)
    mkdir(directory)
  try: 
     results = functional(crystal, outdir=directory, comm=default_comm)
     assert results.success
     # check that final crystals are same as final structures.
     for crystal, structure in zip(results.details.crystal.itervalues(),
                                   results.details.structure.itervalues()):
       estruc = crystal.eval()
       assert all(abs(estruc.cell - structure.cell) < 1e-8)
       assert all(abs( array([a.pos for a in estruc])
                       - array([a.pos for a in structure]) ) < 1e-8)
       assert all( array([a.type for a in estruc])
                   == array([a.type for a in structure]) )
     instruc = results.details.input_crystal.values()[1:] + [results.crystal]
     outstruc = results.details.structure.values()
     for a, b in zip(instruc, outstruc):
       a = a.eval()
       assert all(abs(a.cell - b.cell) < 1e-8)
       assert all(abs( array([v.pos for v in a])
                       - array([v.pos for v in b]) ) < 1e-8)
       assert all( array([v.type for v in a])
                   == array([v.type for v in b]) )
  finally: 
    if directory != '/tmp/test/': rmtree(directory)
    pass
if __name__ == '__main__':
  test()
