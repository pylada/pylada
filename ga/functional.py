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
__all__ = ['Darwin', 'maximize', 'minimize']
from abc import ABCMeta, abstractmethod
from .extract import Extract as GAExtract

class Functional(object):
  """ Genetic Algorithm capable of stop-and-go. """
  __metaclass__ = ABCMeta
  SHELVECAR = 'SHELVECAR'
  """ Name of the file where the history of a GA run is saved. """
  CALCDIR  = 'calcs'
  """ Name of subdirectory where calculations are stored. """
  Extract = GAExtract
  """ Extraction object for GA functional. """

  @abstractmethod
  def random_individual(self):
    """ Returns new random individual. """
    pass
  @abstractmethod
  def jobinator(self, job, indiv):
    """ Fills jobdictionary item for computing individual.
    
        :param job: 
          :py:class:`~pylada.jobfolder.jobfolder.JobFolder` instance which should be
          filled with the right parameters for computing the intermediate
          fitness of an individual.
        :param indiv:
          The individual for which to perform a calculation.
    """
    pass
  @abstractmethod
  def objective(self, extract, indiv):
    """ Returns fitness of an individual.
    
        :param extract: 
          An extraction object which comes as a result of the computation
          created by :py:meth:`~Functional.jobinator`. It should be used
          to actually evaluate the fitness of the individual.
        :param indiv:
          Individual evaluate.
    """
    pass

  def update_fitness(self, population=None):
    """ Updates fitness of whole population. 
    
        This is used to update the fitness of the whole population when it
        changes depending on history, e.g. for X-GSGO.
    """
    return 

  def checkpoints(self):
    """ Whatever needs printing or checking.

        Called at the end of a generational loop.

        :returns: True if the GA should end.
    """
    return False
  
  def go_to_next_iteration(self):
    """ Returns True when going to next generation.

        Not all calculations need be finished.
        It is the job of this function to catch errors in process execution and
        correct them.
    """ 
    from ..process import Fail
    assert hasattr(self.process, '_comm')
    try: self.process.wait()
    except Fail: self._process_errors() 
    return True
  
  def _process_errors(self):
     """ Takes care of errors. """
     offmap = {}
     for i, indiv in enumerate(self.offspring): offmap[str(indiv.uuid)] = i
     uuids = set()
     for name in self.process.errors.iterkeys():
       self.errors.add(name)
       uuid = name.split('/')[-1]
       uuids.add( self.offspring[offmap[uuid]].uuid )
       del self._jobfolder.root[name] 
       print "Calculation {0} failed.".format(name)
     self.process.errors = {}
     self.offspring = [u for u in self.offspring if u.uuid not in uuids]

  @abstractmethod
  def __init__(self, popsize=100, rate=0.1):
    """" Initializes Functional. """
    from .standard import Mating
    super(Functional, self).__init__()

    self.population = []
    """ List of individuals in the population. """
    self.offspring = []
    """ List of uncomputed individuals. """
    self.errors = set()
    """ List of individuals for which calculations failed. """
    self.generation = 0
    """ Current generation.

        The counter is incremented after each successful call to :py:meth:`next`.
    """
    self.matingops = Mating(sequential=False)
    """ Mating operators. """
    self.target_popsize = popsize
    """ Target population size.
    
        This is the size of the population if all goes well.
        Since some calculations may take longuer than others, this size may not
        be attained within the first generation. The method
        :py:meth:`Functional.go_to_next_iteration` decides when to give up on a
        job for now. The job can be incorporated at the next iteration, if it
        finishes in time.
    """
    self.rate = rate
    """ Target turnover rate at each generation. """


  def _restart(self): 
    """ Updates self from SHELVECAR.
    
        Reads the last functional in the shelvecar and updates offspring,
        population, errors, and generation. Nothing else is changed.

        If the shelvecar does not exist, then creates one with current values.
    """
    from os.path import exists
    from shelve import open as shelve_open
    from ..error import internal
    from ..misc import LockFile
    if '_directory' not in self.__dict__:
      raise internal("Directory has not yet been set.")
    with LockFile(self.shelvepath, timeout=10) as lock:
      if not exists(self.shelvepath):
        try:
          shelvecar = shelve_open(self.shelvepath)
          shelvecar['individuals'] = {}
          shelvecar['functionals'] = [self]
          shelvecar['added']       = []
          shelvecar['removed']     = []
          shelvecar['new']         = []
        finally: shelvecar.close()
      else:
        try:
          shelvecar = shelve_open(self.shelvepath, writeback=True)
          previous = shelvecar['functionals']
        finally: shelvecar.close()
        self.offspring  = previous[-1].offspring
        self.population = previous[-1].population
        self.errors     = previous[-1].errors
        self.generation = previous[-1].generation

  def selection(self, dummy, size=2):
    """ Selects individual from the population through deterministic tournament. """
    import random
    list_ = range(len(self.population))
    random.shuffle(list_)
    list_ = list_[:size]
    result = list_[0]
    for b in list_[1:]:
      if self.population[b].fitness < self.population[result].fitness: result = b;
    return result

  @property
  def shelvepath(self):
    """ Directory where output should be found. """
    from os.path import join
    return join(self._directory, self.SHELVECAR)

  @property
  def popsize(self): 
    """ Size of current population. """
    return len(self.population)
  @property
  def nboffspring(self): 
    """ Size of current population. """
    return len(self.offspring)

  def _add_indiv(self, individual, mechanics='mating'):
    """ Adds new individual if it is not taboo. 
    
        Does not update the SHELVECAR file.
    """
    from uuid import uuid4
    if self.is_taboo(individual): return False

    individual.uuid = uuid4()
    individual.conception = self.generation
    individual.mechanics  = mechanics
    self.offspring.append(individual)
    print '{0}: {1}'.format(mechanics, individual)
    return True

  def add_indiv(self, individual, mechanics='manual'):
    """ Adds new individual if it is not taboo. 

        Updates SHELVECAR file. 
    """ 
    if not self._add_indiv(individual, mechanics): return False
    self.update_func()
    return True

  def is_taboo(self, individual):
    """ Returns True if individual is taboo. 


        Generic function checks that the individual is not in offspring or in
        population.
    """
    return individual in self.offspring or individual in self.population

  def mating(self, n=1, itermax=10):
    """ Returns new individual from population. """
    from ..error import internal
    # sanity check
    if len(self.population) == 0: raise internal("No population yet..") 

    dummy = list(self.offspring)
    try:
      # creates new individuals
      for i in xrange(n):
        success = False
        for j in xrange(itermax):
          if self._add_indiv( self.matingops(self) ): 
            success = True
            break
        if not success: raise internal('Could not create new individual.')
    except: # revert changes on failure.
      self.offspring = dummy
      raise

  def new_random_individuals(self, n=1, itermax=10, **kwargs):
    """ Creates a new individual and adds it to the offspring population. """
    # creates new individuals
    dummy = list(self.offspring)
    try:
      for i in xrange(n):
        success = False
        for j in xrange(itermax):
          if self._add_indiv( self.random_individual(**kwargs), mechanics='random' ): 
            success = True
            break
        if not success: raise RuntimeError('Could not create new individual.')
    except: # revert changes on failure.
      self.offspring = dummy
      raise

  @property
  def calcdir(self):
    """ Returns current jobdict path. """
    from os.path import join
    return join(self._directory, self.CALCDIR)

  @property
  def jobfolder(self):
    """ Returns current jobfolder.
    
        The current job-folder holds all the currently running jobs.
        It is a property because it is not saved to disk when pickle. 
        It should be created anew each time GA runs.
    """
    from ..jobfolder import JobFolder
    if '_jobfolder' not in self.__dict__:
      self._jobfolder = JobFolder()
      for indiv in self.offspring: 
        job = self._jobfolder / str(indiv.conception) / str(indiv.uuid)
        job.uuid = indiv.uuid
        job._indiv = indiv
        self.jobinator(job, indiv)
    return self._jobfolder

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
    from ..process import JobFolderProcess
    if '_process' not in self.__dict__: 
      self._process = JobFolderProcess(self.jobfolder, self.calcdir, keepalive=True)
    return self._process

  def __call__(self, outdir, comm): 
    """ Performs GA run until exhaustion sets in. """
    import sys
    from itertools import chain
    from os.path import join
    from shelve import open as shelve_open
    from ..misc import RelativePath, LockFile

    # sets ouput directory for this run.
    self._directory = RelativePath(outdir).path
    # reloads population, offspring, errors, generation from previous run
    # this has to come after setting the directory.
    self._restart()

    # try loops resets old stdout, stderr
    oldouterr = sys.stdout, sys.stderr
    try: 
      sys.stdout = open(join(self._directory, "out"), 'a')
      sys.stderr = open(join(self._directory, "err"), 'a')
     
      print 'COMMUNICATOR', comm 

      # first creates a random population if necessary.
      if len(self.population) == 0 and len(self.offspring) == 0: 
        self.new_random_individuals(n=self.target_popsize)
        self.update_func()
        sys.stdout.flush()

      # then  creates the jobfolder and start the process.
      self.process.start(comm)
      assert hasattr(self.process, '_comm')

      # finally, loop till eternity comes.
      while True:
        if not self.go_to_next_iteration(): continue

        # Now we create new generation. 
        # First we get those individuals that were successfuly computed.
        successfuls = self.find_successfuls()

        # update fitness of whole population, in case it depends upon whole population.
        iterator = (self.offspring[i] for i in successfuls.itervalues())
        self.update_fitness(chain(iterator, self.population))

        # prints out successful individuals.
        if len(successfuls) != 0:
          print "Newly calculated individuals. """
          for index in successfuls.itervalues():
            print " {0.uuid}(conception {0.conception}): {0.fitness}"\
                  .format(self.offspring[index])
          sys.stdout.flush()
         
        # removes worst individuals.
        removed, population, offspring = self.survival(successfuls)
        self.population = population
        self.offspring = offspring

        # increment generation count.
        self.generation += 1

        # create new offspring via mating
        nboffspring = int(self.target_popsize*self.rate+0.1) - len(offspring) \
                      + (self.target_popsize - self.popsize)
        self.mating(nboffspring)

        # udpate jobfolder
        self._update_jobfolder(removed, offspring, population)

        # update shelvecar
        with LockFile(self.shelvepath, timeout=10) as lock:
          try:
            shelvecar = shelve_open(self.shelvepath, writeback=True)
            individuals = shelvecar['individuals']
            for individual in self.population: 
              if not hasattr(individual, 'birth'):
                individual.birth = self.generation
              individuals[str(individual.uuid)] = individual
            for uuid in removed: individuals[str(uuid)].removed = self.generation
            shelvecar['individuals'] = individuals
            shelvecar['functionals'].append(self)
            shelvecar['removed'].append(removed)
            shelvecar['added'].append(successfuls.keys())
            shelvecar['new'].append( [ u.uuid for u in self.offspring\
                                       if u.conception == self.generation-1] )
          finally: shelvecar.close()
        print "\n Starting generation {0}.".format(self.generation)
        sys.stdout.flush()

        if  self.checkpoints(): break

    finally:
      sys.stdout.close()
      sys.stderr.close()
      sys.stdout, sys.stderr = oldouterr 
      if '_process' in self.__dict__:
        self._process.terminate()
        del self._process
      self.__dict__.pop('_directory', None)
      self.__dict__.pop('_jobfolder', None)
    # returns an extraction object. 
    # Might never get here, but that's ok.
    return self.Extract(outdir)

  def find_successfuls(self):
    """ Mapping of successful calculations to individuals. """
    from os.path import join
    from ..jobfolder import MassExtract
    # creates mass extraction object with fake jobfolder file, but real
    # jobfolder.
    collect = MassExtract(join(self.calcdir, 'folder'))
    collect._jobfolder = self.jobfolder
    # creates map from uuid to index in list of offspring.
    offmap = {}
    for i, indiv in enumerate(self.offspring): offmap[indiv.uuid] = i
    # now check successful runs, gets fitness for those, and create mapping
    # from names to index in offspring list.
    successfuls = {}
    for name, extractor in collect.iteritems(): 
      if extractor.success == False: continue

      # try and compute that individual's fitness
      uuid = self.jobfolder[name].uuid
      indiv = self.offspring[offmap[uuid]]
      indiv.fitness = self.objective(extractor, indiv)
      successfuls[name] = offmap[uuid]
    return successfuls

  def survival(self, successfuls):
    """ Removes worst individuals. """
    from operator import attrgetter
    # reates list of offspring
    offspring = [ indiv for indiv in self.offspring \
                  if "/{0.conception}/{0.uuid}/".format(indiv) \
                     not in successfuls ]
    if len(successfuls) == 0: return [], list(self.population), offspring
    # Removes from population the worst individuals.
    howmany = 0 if self.popsize == 0 \
              else len(successfuls) - (self.target_popsize - self.popsize)
    if howmany < 0: howmany = 0

    # Minimizing: remove individuals with lowest fitness.
    population = sorted(self.population, key=attrgetter('fitness'))
    removed = [u.uuid for u in population[len(population)-len(successfuls):]]
    # and creates new population
    population = population[:len(population) - howmany] \
                 + [ self.offspring[i] for i in successfuls.itervalues() ]
    return removed, population, offspring


  def update_func(self):
    """ Updates functional in SHELVECAR file. """
    from shelve import open as shelve_open
    from ..misc import LockFile
    with LockFile(self.shelvepath, timeout=10) as file:
      try:
        shelvecar = shelve_open(self.shelvepath, writeback=True)
        shelvecar['functionals'][-1] = self
      finally: shelvecar.close()

  def _update_jobfolder(self, removed, successfuls, offspring):
     """ Updates jobfolder and process. """
     del self._jobfolder
     self.process.update(self.jobfolder, deleteold=True)

  def __getstate__(self):
    """ Returns current state.

        Current state does not include directory, jobfolder, or process.
    """ 
    result = self.__dict__.copy()
    result.pop('_process', None)
    result.pop('_jobfolder', None)
    result.pop('_directory', None)
    return result

  def __setstate__(self, value):
    """ Sets state. """
    # removes current process, directory, jobfolder.
    if '_process' in self.__dict__:
      self._process.terminate()
      del self._process
    self.__dict__.pop('_directory', None)
    self.__dict__.pop('_jobfolder', None)
    self.__dict__.update(value)
