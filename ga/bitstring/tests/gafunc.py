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

from pylada.ga.functional import Functional as GAFunctional
class Functional(GAFunctional):
  """ GA bitstring functional. """
  def __init__(self, program, size=10, maxgen=-1, **kwargs):
    """ Initializes functional. """
    from pylada.ga import bitstring
    super(Functional, self).__init__(**kwargs)
    self.size = size
    self.program = program
    self.matingops.add( bitstring.Crossover(rate=0.25), rate=0.8 )
    self.matingops.add( bitstring.Mutation(rate=3e0/float(bitstring.Individual.size)), rate=0.2 )
    self.maxgen = maxgen

  def random_individual(self):
    """ Returns a new random individual. """
    from pylada.ga.bitstring import Individual
    return Individual(size=self.size)

  def jobinator(self, job, indiv):
    """ Initializes a new job from an individual. """
    from functional import SerialFunctional
    from random import random
    job.functional = SerialFunctional(self.program, order=sum(indiv.genes))
    if random() > 0.8: # generate artificial failure
      job.functional.order = 666
  def objective(self, extract, indiv):
    """ Returns error. """
    return extract.error

  def checkpoints(self):
    """ Returns true when result is in population. """
    from numpy import all
    if self.generation > self.maxgen and self.maxgen >= 0: return True
    for indiv in self.population:
      if all(indiv.genes == 1): return True
    return False
  def is_taboo(self, individual, mechanics='mating'): return False

  @property 
  def process(self):
    """ Current process. 

        The process is an object which manages launching actual calculation.
        It basically takes a job-folder and launches calculations as
        appropriate. 

        By default, it creates a
        :py:class:`~pylada.process.jobfolder.JobFolderProcess`. It might be a
        good idea to use something more specific, such as a
        :py:class:`~pylada.process.pool.PoolProcess` instantiated with the
        correct ressource allocation function.

        The process is not saved when self is pickled. It is created anew each
        time this functional runs.
    """ 
    from pylada.process import JobFolderProcess
    if '_process' not in self.__dict__: 
      self._process = JobFolderProcess(self.jobfolder, self.calcdir, nbpools=4, keepalive=True)
    return self._process
