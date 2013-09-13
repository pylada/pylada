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

def create(input='input.py'):
  """ Load ga into ipython. """
  from os.path import exists
  from pickle import dump
  from pylada.vasp import read_input
  from pylada.ga.xgsgo.functional import Functional
  from pylada.jobfolder import JobFolder

  input = read_input(input)


  root = JobFolder()
  for trial in xrange(input.trials):
    folder = root / "results" / str(trial)
    folder.functional = Functional( functional        = input.vasp, 
                                    species           = input.species,
                                    natoms            = input.natoms,
                                    rate              = input.rate,
                                    popsize           = input.popsize,
                                    cns_rate          = input.cns_rate,
                                    mix_atoms_rate    = input.mix_atoms_rate,
                                    mix_poscar_rate   = input.mix_poscar_rate,
                                    mutation_rate = input.mutation_rate )
  return root
