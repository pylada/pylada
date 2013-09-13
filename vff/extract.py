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

""" Subpackage containing extraction methods for CRYSTAL output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']

from ..tools import make_cached
from ..tools.extract import search_factory, AbstractExtractBase
from ..error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'vff.out')

class ExtractBase(object):
  """ Implementation class for extracting data from VFF output """

  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()


  @property
  def energy(self):
    """ Strain energy. """
    return self.structure.energy
  @property
  def stress(self):
    """ Stress tensor. """
    return self.structure.stress
  @property
  def forces(self):
    """ Forces a an nx3 array. """
    from numpy import array
    return array([u.gradient for u in self.structure])                         \
           * self.structure[0].gradient.units
  @property
  @make_cached
  def start_date(self):
    """ Starting date of the calculation. """
    from datetime import datetime
    regex = self._find_first_STDOUT('^Start date:\s*(.*)$')
    if regex is None: raise GrepError('Could not find starting date.')
    regex = regex.group(1).split()
    regex = list(regex[0].split('-')) + list(regex[1].split(':'))
    regex = regex[:-1] + list(regex[-1].split('.'))
    regex = [int(u) for u in regex]
    return datetime(*regex)
  @property
  @make_cached
  def end_date(self):
    """ Starting date of the calculation. """
    from datetime import datetime
    regex = self._find_first_STDOUT('^End date:\s*(.*)$')
    if regex is None: raise GrepError('Could not find end date.')
    regex = regex.group(1).split()
    regex = list(regex[0].split('-')) + list(regex[1].split(':'))
    regex = regex[:-1] + list(regex[-1].split('.'))
    regex = [int(u) for u in regex]
    return datetime(*regex)

  @property
  @make_cached
  def input_structure(self):
    """ Structure prior to optimization. """
    from . import exec_input
    script = self._extract_script("INITIAL STRUCTURE")
    return exec_input(script).structure

  @property
  @make_cached
  def structure(self):
    """ Structure after optimization. """
    from . import exec_input
    script = self._extract_script("STRUCTURE")
    return exec_input(script).structure

  @property
  @make_cached
  def functional(self):
    """ Functional used for optimization. """
    from . import exec_input
    script = self._extract_script("FUNCTIONAL")
    return exec_input(script).functional

  @property
  @make_cached
  def optimize(self):
    """ Optimization output. """
    from collections import namedtuple
    from re import compile
    import numpy
    script = self._extract_script("MINIMIZATION")
    if script is None:
      raise GrepError('No information about optimization.')
    script = script.splitlines()
    allinputs = {}
    current = None
    pattern = compile("^\s*(\S+)\s*:")
    for line in script:
      regex = pattern.search(line)
      if regex is not None:
        current = regex.group(1)
        line = line[line.find(':')+1:]
        allinputs[current] = line
      elif current is not None:
        allinputs[current] += '\n' + line
    Result = namedtuple('VffOptimization', allinputs.keys())
    global_dict = numpy.__dict__.copy()
    for key in allinputs.iterkeys():
      try: dummy = eval(allinputs[key], global_dict)
      except: pass
      else: allinputs[key] = dummy
    return Result(**allinputs)

  def _extract_script(self, pattern):
    """ Returns string between two '#+' patterns. """
    from re import compile
    with self.__stdout__() as file:
      regex = compile('^#+ {0} #+$'.format(pattern))
      found = False
      for line in file:
        if regex.search(line) != None: found = True; break
      if not found: return None
      result = ""
      regex = compile('^#+ END {0} #+$'.format(pattern))
      found = False
      for line in file:
        if regex.search(line) != None: found = True; break
        result += line
      return result if found else None
     



class Extract(AbstractExtractBase, OutputSearchMixin, ExtractBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object. 
    
        :param directory: 
          Directory where the OUTCAR resides. 
          It may also be the path to an OUTCAR itself, if the file is not
          actually called OUTCAR.
    """
    from os.path import exists, isdir, basename, dirname
    from pylada.misc import RelativePath
       
    self.STDOUT = 'vff.out'
    """ Name of file to grep. """
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        self.STDOUT = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractBase.__init__(self)
    OutputSearchMixin.__init__(self)

  @property
  def success(self):
    try: return self.optimize.success
    except: return False

  def iterfiles(self, **kwargs):
    """ iterates over input/output files. 
    
        :param bool input: Include INCAR file
        :param bool wavefunctions: Include WAVECAR file
        :param bool structure: Include POSCAR file
    """
    from os.path import exists, join
    file = join(self.directory, self.STDOUT)
    if exists(file): yield file

del make_cached
del search_factory
