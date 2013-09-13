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

species = ['C', 'N']
""" Species involved in the optimization. """
anions  = ['N']
""" Anions, if any, among the species.

    All anions must be included in the species.
"""
natoms  = (2, 8)
""" Min, Max number of atoms. """
rate    = 0.5
""" Offspring mutation rate.  """
popsize = 10
""" Target population size. """
cns_rate = 0.8
""" Cut and splice crossover rate. """
mutation_rate = 0.1
""" Cut and splice crossover rate. """
mix_atoms_rate = -1
""" mix-atom crossover rate. """
mix_poscar_rate = -1
""" mix-poscar crossover rate. """
trials = 4
""" Number of independent trials. """

first_trial = {'kpoints': '\n0\nAuto\n1\n'}
vasp = Relax(encut=1.0, kpoints='\n0\nAuto\n10', relaxation='volume ionic cellshape')
vasp.add_specie = 'C', 'pseudos/C'
vasp.add_specie = 'N', 'pseudos/N'
