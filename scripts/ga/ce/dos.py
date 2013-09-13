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

""" Creates job-dictionary of GA jobs. """
def gajobs(path, inputpath = "input.py", collector=None):
  from weakref import proxy

  import IPython.ipapi
  from numpy import array, arange, any, isnan
  from quantities import eV

  from pylada.ga.functional import maximize, minimize
  from pylada.ga.ce.functional import Darwin as Functional
  from pylada.ga.ce.evaluator  import LocalSearchEvaluator as Evaluator
  from pylada.jobs import JobFolder
  from pylada.escan import MassExtract
  from pylada.ce import read_input

  # reads input file.
  input = read_input(inputpath, {"maximize": maximize, "minimize": minimize, "eV": eV, "arange": arange})

  if collector is None: collector = MassExtract(input.directory)
  extractor  = [extract for extract in collector.values() if extract.success]
  structures = [extract.input_structure for extract in extractor]
  for e in extractor: assert not any(isnan(e.eigenvalues))
  dos_values = array([extract.dos(input.energies, input.sigma) for extract in extractor]).T
  assert not any(isnan(dos_values))
  evaluator  = Evaluator(input.lattice, structures, energies=dos_values[0], **input.clusters)


  jobfolder = JobFolder()
  for energy, energies in zip(input.energies, dos_values):
    for trial in input.trials:
      kwargs = { "popsize": input.population_size, "rate": input.offspring_rate,
                 "max_gen": input.max_generations, "crossover_rate": input.crossover_rate, 
                 "history": input.history, "histlimit": input.histlimit,
                 "comparison": input.comparison, "rootworkdir": "$SCRATCH",
                 "alwayson": input.alwayson, "alwaysoff": input.alwaysoff }
      functional = Functional(evaluator, **kwargs)
      evaluator.darwin = proxy(functional)

      gajob = jobfolder / "dos_{0}/trial_{1}".format(energy, trial)
      print "dos_{0}/trial_{1}".format(energy, trial)
      gajob.functional = functional 
      gajob.jobparams["energies"] = energies

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobfolder"] = jobfolder
  ip.magic("savejobs " + path)
