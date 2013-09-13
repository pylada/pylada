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

""" GA as a functional, for use with pyladajobs. """
__docformat__ = "restructuredtext en"
__all__ = ['Darwin']
from ..functional import Darwin as StandardDarwin

class Darwin(StandardDarwin):
  """ Main GA functional for CE. """
  def __init__(self, evaluator, **kwargs):
    """ Initializes a GA functional. 
         
        :Parameters:
          evaluator : `pylada.ga.escan.elemental.Bandgap`
            Functional which uses vff and escan for evaluating physical properties.
        :Kwarg rate: 
          Offspring rate. Defaults to ``max(0.1, 1/popsize)``.
        :Kwarg popsize: 
          Population size. Defaults to 100.
        :Kwarg maxgen: 
          Maximum number of generations. Defaults to 0.
        :Kwarg current_gen:
          Current generation. Defaults to 0 or to whatever is in the restart.
        :Kwarg alwayson:
          Always leave these figures on.
        :Kwarg alwaysoff:
          Always leave these figures off.
    """
    from . import Crossover, Mutation
    from ..standard import Mating
    # remove pools.
    kwargs.pop("pools", None)
    kwargs.pop("rootworkdir", None)
    alwayson = kwargs.pop("alwayson", None)
    alwaysoff = kwargs.pop("alwaysoff", None)

    super(Darwin, self).__init__(evaluator, **kwargs)

    # creates mating operation.
    self.crossover_rate = kwargs.pop('crossover_rate', 0.8)
    """ Rate for crossover over other operations. """
    self.swap_rate = kwargs.pop('swap_rate', 0.8)
    """ Rate for swap-like mutationns over other operations. """
    self.matingops = Mating(sequential=False)
    self.matingops.add(Crossover(alwayson=alwayson, alwaysoff=alwaysoff), self.crossover_rate)
    self.matingops.add( Mutation(len(evaluator.clusters), alwayson=alwayson, alwaysoff=alwaysoff),\
                        1-self.crossover_rate )
    self._alwayson, self._alwaysoff = set(), set()
    self.alwayson = alwayson
    """ Clusters which are always on. """
    self.alwaysoff = alwaysoff
    """ Clusters which are always off. """

  @property
  def alwayson(self):
    """ Clusters which are always on. """
    return self._alwayson
  @property
  def alwaysoff(self):
    """ Clusters which are always off. """
    return self._alwaysoff
  @alwayson.setter
  def alwayson(self, value):
    self._alwayson = set(value) if value is not None else set()
    for function, dummy, dummy in self.matingops.operators:
      if hasattr(function, "alwayson"): setattr(function, "alwayson", self._alwayson)
    self.evaluator.exclude = self._alwayson | self._alwaysoff
  @alwaysoff.setter
  def alwaysoff(self, value):
    self._alwaysoff = set(value) if value is not None else set()
    for function, dummy, dummy in self.matingops.operators:
      if hasattr(function, "alwaysoff"): setattr(function, "alwaysoff", self._alwaysoff)
    self.evaluator.exclude = self._alwayson | self._alwaysoff


  def taboo(self, indiv):
    """ Makes sure that sets are not empty. """
    if len(indiv.genes) == 0: return True
    return super(Darwin, self).taboo(indiv)

  def mating(self):
    """ Calls mating operations. """
    return self.matingops(self)

  def Individual(self):
    """ Generates an individual (with mpi) """
    from . import Individual
    N = len(self.evaluator.clusters)
    result = Individual(N, mean=int(0.5 * N), alwayson=self.alwayson, alwaysoff=self.alwaysoff)
    return result

  def _further_setup(self):
    """ Last minute setup. """
    from weakref import proxy
    self.pools = self.comm.size
    self.evaluator.darwin = proxy(self)


  def evaluation(self):
    """ Evaluates population. """
    from ..standard import mpi_population_evaluation
    mpi_population_evaluation(self, self.evaluator, self.pools, self.comm)

  @property
  def rootworkdir(self):
    """ Root of the working directory where actual calculations are carried out. """
    raise RuntimeError("No rootworkdir necessary for CE.")

  @rootworkdir.setter
  def rootworkdir(self, value): pass
