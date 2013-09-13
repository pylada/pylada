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

""" Extraction module for escan elemental GA. """
__docformat__ = "restructuredtext en"
__all__  = ['Extract']

class Extract(object):
  """ Extraction class for GA Functional. """
  def __init__(self, directory=None):
    """ Initializes an extraction object. 
  
        :param directory: 
          Directory where the SHELVECAR can be found.
          Defaults to current working directory.

        :raises RuntimeError: if no appropriate SHELVECAR can be found.
    """
    from os.path import exists, isfile, isdir
    from shelve import open as shelve_open
    from ..misc import RelativePath, LockFile
    super(Extract, self).__init__()

    self._directory = RelativePath(directory)
    if not exists(self.directory):return
    if not isdir(self.directory): return
    if not exists(self.shelvepath): return
    if not isfile(self.shelvepath): return

    with LockFile(self.shelvepath, timeout=10) as lock:
      shelve = shelve_open(self.shelvepath)
      try:
        if set(shelve.keys()) !=  set(['individuals', 'functionals', 'removed', 'added', 'new']):
          raise RuntimeError('{0} is not a GA SHELVECAR file.'.format(self.shelvepath))
      finally: shelve.close()
  @property
  def success(self): return False

  @property
  def directory(self):
    """ Directory where output should be found. """
    return self._directory.path
  @property
  def shelvepath(self):
    """ Directory where output should be found. """
    from os.path import join
    from .functional import Functional
    return join(self.directory, Functional.SHELVECAR)

  @property
  def functionals(self):
    """ Functionals at each generation. """
    from shelve import open as shelve_open
    from ..misc import LockFile
    with LockFile(self.shelvepath, timeout=10) as lock:
      try:
        shelve = shelve_open(self.shelvepath)
        return shelve['functionals']
      finally: shelve.close()
  @property 
  def populations(self):
    """ Populations at each generation. """
    return [func.population for func in self.functionals]

  @property
  def functional(self):
    """ Return current functional. """
    return self.functionals[-1]

  @property
  def generation(self):
    """ Returns current generation. """
    return self.functional.generation
  @property
  def population(self):
    """ Returns current population. """
    return self.functional.population
  @property
  def offspring(self):
    """ Returns current population. """
    return self.functional.offspring
    

  @property 
  def added(self):
    """ Individuals added at each generation. """
    from shelve import open as shelve_open
    from ..misc import LockFile
    with LockFile(self.shelvepath, timeout=10) as lock:
      try:
        shelve = shelve_open(self.shelvepath)
        result = []
        for add in shelve['added']:
          result.extend(shelve['individuals'][i] for i in add)
        return result
      finally: shelve.close()

  @property 
  def removed(self):
    """ Individuals removed at each generation. """
    from shelve import open as shelve_open
    from ..misc import LockFile
    with LockFile(self.shelvepath, timeout=10) as lock:
      try:
        shelve = shelve_open(self.shelvepath)
        result = []
        for r in shelve['removed']:
          result.extend(shelve['individuals'][i] for i in r)
        return result
      finally: shelve.close()

  def best(self, n=10): 
    """ Current n best individuals. 

        :param int n: 
           Number of individuals to print.
    """
    from operator import itemgetter
    from ..jobfolder.forwarding_dict import ForwardingDict

    individuals = self.individuals

    fitness = sorted( [(name, u.fitness) for name, u in individuals.iteritems()],\
                      key=itemgetter(1) )[:n]
    result = ForwardingDict(readonly=True, ordered=False)
    for key, value in fitness: result[key] = individuals[key]
    return result
    
  @property
  def individuals(self):
    """ Dictionary containing all individuals. """
    from shelve import open as shelve_open
    from ..misc import LockFile
    with LockFile(self.shelvepath, timeout=10) as lock:
      try: shelve = shelve_open(self.shelvepath)
      except: raise
      else: return shelve['individuals']
      finally: shelve.close()


