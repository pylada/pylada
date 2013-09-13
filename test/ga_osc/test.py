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
def gajobs(path, inputpath = "input.py"):
  from copy import deepcopy
  import IPython.ipapi

  from pylada.ga.escan.elemental.functional import Darwin as Functional
  from pylada.ga.escan.elemental.evaluator  import Dipole as Evaluator
  from pylada.ga.escan.elemental            import LayeredConverter as Converter
  from pylada.jobs import JobFolder
  from pylada.escan import read_input

  # reads input file.
  input = read_input(inputpath)

  jobfolder = JobFolder()
  for scale in input.scales:
    for p in input.periods: 
      for trial in input.trials:
        supercell = input.supercell(p)

        escan = deepcopy(input.escan)
        escan.fft_mesh = input.fftmesh(supercell)
        escan.vff.lattice.scale = scale

        converter = Converter(supercell, lattice=escan.vff.lattice)
        
        evaluator = Evaluator(converter = converter, escan = escan, **input.kwargs)
        evaluator._outdir.envvar = "$SCRATCH"

        functional = Functional(evaluator)
        functional.popsize     = input.population_size
        functional.max_gen     = input.max_generations
        functional.rate        = input.offspring_rate
        functional.mean_conc   = input.mean_conc
        functional.stddev_conc = input.stddev_conc
        functional.pools       = input.pools

        gajob = jobfolder / "scale_{0:.2f}/period_{1}/trial_{2}".format(scale, p, trial)
        gajob.functional = functional 

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobfolder"] = jobfolder
  ip.magic("savejobs " + path)
