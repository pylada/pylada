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

from ..functional import Functional as GAFunctional
from .objective import objective as objective_xgsgo
class Functional(GAFunctional):
  """ XGSGO GA functional. """
  def __init__(self, functional, species, natoms=None, anions=None,
               cns_rate=1e0, mix_atoms_rate=1e0, mix_poscar_rate=1e0,
               mutation_rate=0.5, **kwargs): 
    from .operators import cut_and_splice, mix_atoms, mix_poscars,\
                           jiggle_structure
    super(Functional, self).__init__()
    if cns_rate > 0: self.matingops.add(cut_and_splice, rate=cns_rate)
    if mix_atoms_rate > 0: self.matingops.add(mix_atoms, rate=mix_atoms_rate)
    if mix_poscar_rate > 0:
      self.matingops.add(mix_poscars, rate=mix_poscar_rate)
    if mutation_rate > 0:
      self.matingops.add(jiggle_structure, rate=mutation_rate)

    self.functional = functional
    """ Functional with which to get total energy. """
    self.species    = list(species)
    """ Species in GA optimization. """
    self.anions     = anions
    """ Anionic species in the calculation. """
    self.natoms     = (2, 20) if natoms is None else (min(natoms), max(natoms))
    """ Maximum number of atoms. """
    if self.natoms[0] < 1: self.natoms = (1, self.natoms[1])

  def random_individual(self):
    """ Returns a new random individual. """
    from random import randint, shuffle
    from .initialization import random_structure, populate, populate_anion_cation
    N = randint(*self.natoms)
    result = random_structure(N)
    natoms = []
    for i in xrange(len(self.species)):
      natoms.append(N if N < 2 else randint(1, N))
      N -= natoms[-1]
    if N != 0: natoms[randint(0, len(natoms)-1)] += N
    shuffle(natoms)
    species = {}
    for s, n in zip(self.species, natoms): species[s] = n
    if self.anions is not None and len(self.anions) > 0:
      populate_anion_cation(result, species, self.anions)
    else: populate(result, species)
    return result

  def jobinator(self, job, indiv):
    """ Initializes a new job from an individual. """
    job.functional = self.functional
    job.params['structure'] = indiv
      
  objective = objective_xgsgo
  """ Stores in individual information needed for convex-hull. 

      Does not actuall compute fitness, since we need all individuals for that.
  """
  def is_taboo(self, individual):
    """ Checks whether an individual is taboo. """
    from .initialization import taboo
    if len(individual) < self.natoms[0] or len(individual) > self.natoms[1]: 
      return False
    return taboo(individual, same_first_neigh=-1)

  @property 
  def process(self):
    """ Current process. 

        The process is an object which manages launching actual calculation.
        It basically takes a job-folder and launches calculations as
        appropriate. 

        By default, it creates a
        :py:class:`~pylada.process.pool.PoolProcess` where each job is allocated
        the N procs, where N is the even number closest from below to the
        number of atoms in the stucturel.

        The process is not saved when self is pickled. It is created anew each
        time this functional runs.
    """ 
    from pylada.process import PoolProcess
    def nbprocs(job):
      return len(job.structure) - len(job.structure) % 2
    if '_process' not in self.__dict__: 
      self._process = PoolProcess(self.jobfolder, self.calcdir, nbprocs, keepalive=True)
    return self._process

  def __call__(self, *args, **kwargs):
    """ Performs GA. """
    self.nbprocs = kwargs['comm']['n']
    return super(Functional, self).__call__(*args, **kwargs)

  def go_to_next_iteration(self):
    """ Returns True when going to next generation.

        Not all calculations need be finished.
        It is the job of this function to catch errors in process execution and
        correct them.
    """ 
    from time import sleep
    from ...process import Fail
    
    fewjobsleft = max(len(self.offspring) * 0.1, 0)
    lotsofcpusidle = min(self.process._alloc.values())
    while self.process.nbjobsleft > fewjobsleft \
          and self.process._comm['n'] < lotsofcpusidle:
      try:
        if self.process.poll(): break
      except Fail: break
      else: sleep(5)
    self._process_errors()
    return True
