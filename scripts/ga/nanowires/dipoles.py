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

def create_jobs(path, inputpath="input.py", **kwargs):
  """ Creates GA-jobs. 
  
      :Parameters:
        path 
          Path where the job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. Calculations will be
          performed in same directory as this file.
        inputpath
          Path to an input file. Defaults to input.py. 
        kwargs
          Any keyword/value pair to take precedence over anything in the input
          file.
  """
  from IPython.ipapi import get as get_ipy

  from pylada.jobs import JobFolder
  from pylada.ga.escan.nanowires import read_input
  from pylada.ga.escan.nanowires.converter import Converter
  from pylada.ga.escan.nanowires.functional import Darwin

  input = read_input(inputpath)

  lattice = input.vff.lattice
  types = set([u for a in lattice.sites for u in a.type]) - set([input.passivant])

  jobfolder = JobFolder()
  for core_type in input.core_types:
    for nmin, nmax in input.ranges:
      for trial in xrange(input.nb_trials):
        # create conversion object from bitstring to nanowires and back.
        converter = Converter( input.vff.lattice, growth=input.growth_direction,
                               core_radius=input.core_radius, core_type=core_type,
                               types = list(types), thickness=input.thickness,
                               sep = input.__dict__.get('separation', 1) )
        # create objective fucntion.
        evaluator = input.Evaluator( converter, input.escan,
                                     degeneracy = input.__dict__.get('degeneracy', 1e-3) )
        
        # create GA functional itself.
        functional = Darwin(evaluator, nmin=nmin, nmax=nmax, **input.__dict__)
        functional.popsize     = input.population_size
        functional.max_gen     = input.max_generations
        functional.rate        = input.offspring_rate
        functional.pools       = input.pools
        if input.rootworkdir is not None: functional.rootworkdir = input.rootworkdir
    
        gajob = jobfolder / "{0}_core/{2}_{3}/trial_{1}".format(core_type, trial, nmin, nmax)
        gajob.functional = functional



  # saves jobs.
  ip = get_ipy()
  ip.user_ns["current_jobfolder"] = jobfolder
  ip.magic("savejobs " + path)

