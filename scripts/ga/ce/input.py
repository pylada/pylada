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

""" Input script to CE-dos GA. """
# vff parameters
from pylada.crystal.binary import zinc_blende
lattice = zinc_blende()
for site in lattice.sites: site.type = "Si", "Ge"

directory = "../start"
""" Input directories with calculations. """

clusters = {"J0":True, "J1":True, "B2":12, "B3":4, "B4":2, "B5":2}

# GA parameters.
population_size = 8
""" Size of the GA population """
max_generations = -1
""" Maximum number of generations """
offspring_rate  = 0.25
""" Rate at which offspring are created. """
crossover_rate = 0.75
""" Rate of crossover over other operations. """
trials = range(3)
""" Number of trial GA runs. """
comparison = minimize
""" What we are trying to do with the fitness. """
history = True
""" Whether to keep track of candidates in order never to visit twice. """
histlimit = int(1.5*population_size)
""" Maximum number of individuals to keep in history. """

energies = arange(11, dtype="float64") / 10. * 5 - 2.5
""" Energies over which to perform DOS. """ 
sigma = 0.1 * eV
""" Gaussian smearing for density of states. """

maxiter = -1
""" Maximum of iterations during local search. """
maxeval = 50
""" Maximum of evaluations during local search. """

alwayson = [0, 1]
""" Always leave these figures on. """
alwaysoff = []
""" Always leave these figures off. """
