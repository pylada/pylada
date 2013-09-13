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

""" Sub-package for gulp extraction object. """
__docformat__ = "restructuredtext en"
__all__ = ['Extract']
from ..tools import make_cached
from ..tools.extract import search_factory, AbstractExtractBase
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'gulp.gout')

class ExtractBase(object):
  """ Greps some stuff from GULP output file. """
  def __init__(self):
    """ Initializes extraction class. """
    super(ExtractBase, self).__init__()

  @property
  @make_cached
  def optimize(self):
    """ True if optimizing structure. """
    return self._find_first_STDOUT('\*\s*optimise\s*\-\s*perform') is not None
  @property
  @make_cached
  def conp(self):
    """ True if optimizing at constant pressure. """
    return self._find_first_STDOUT('\*\s*conp\s*\-\s*constant') is not None
  @property
  @make_cached
  def conv(self):
    """ True if optimizing at constant volume. """
    return self._find_first_STDOUT('\*\s*conv\s*\-\s*constant') is not None
  @property
  @make_cached
  def qeq(self):
    """ True if using QEe. """
    return self._find_first_STDOUT('\*\s*QEq\s*\-\s*electroneg') is not None
  @property
  @make_cached
  def cell(self):
    """ Greps output cell. """
    # static case
    if not self.optimize: return self.input_cell
    if self.conv: return self.input_cell

    from re import M
    from numpy import array
    from quantities import angstrom
    from ..error import GrepError
    expr = '\s*Final\s+Cartesian\s+lattice\s+vectors\s*\(Angstroms\)\s*:\s*'   \
           '\\n\s*\\n'                                                         \
           '(?:\s*(?:(\S*)\s*){3}\s*\\n){3}'
    result = self._find_first_STDOUT(expr, M)
    if result is None: raise GrepError('Could not grep cell')
    return array(result.group(0).split()[-9:], dtype='float64')                \
                .reshape(3,3) * angstrom

  @property
  @make_cached
  def dimensionality(self):
    """ Dimensionality of the structure. """
    from ..error import GrepError
    result = self._find_first_STDOUT('Dimensionality\s*=\s*(0|1|2|3)')
    if result is None: raise GrepError('Could not determine dimensionality.')
    return int(result.group(1))
  @property
  @make_cached
  def input_cell(self):
    """ Greps input cell. """
    from numpy import identity
    if self.dimensionality == 3: return self._input_cell3d()
    elif self.dimensionality == 2: return self._input_cell2d()
    return identity(3, dtype='float64') * 500.0

  def _input_cell3d(self):
    """ Greps input cell from crystal calculation. """
    from re import M
    from numpy import array
    from quantities import angstrom
    from ..error import GrepError
    expr = '\s*Cartesian\s+lattice\s+vectors\s*\(Angstroms\)\s*:\s*'           \
           '\\n\s*\\n'                                                         \
           '(?:\s*(?:(\S*)\s*){3}\s*\\n){3}'
    result = self._find_first_STDOUT(expr, M)
    if result is None: raise GrepError('Could not grep cell')
    return array(result.group(0).split()[-9:], dtype='float64')                \
                .reshape(3,3) * angstrom

  def _input_cell2d(self):
    """ Greps input cell from thin film calculation. """
    from re import M
    from numpy import array, sqrt, cross, sum
    from quantities import angstrom
    from ..error import GrepError
    expr = '\s*Surface\s*Cartesian\s+vectors\s*\(Angstroms\)\s*:\s*'           \
           '\\n\s*\\n'                                                         \
           '(?:\s*(?:(\S*)\s*){3}\s*\\n){2}'
    result = self._find_first_STDOUT(expr, M)
    if result is None: raise GrepError('Could not grep cell')
    a0 = array(result.group(0).split()[-6:-3], dtype='float64')
    a1 = array(result.group(0).split()[-3:], dtype='float64')
    a2 = cross(a0, a1)
    a2 = a2 * 500e0 / sqrt(sum(a2*a2))
    return array([a0, a1, a2]).T * angstrom

  @property 
  @make_cached
  def structure(self):
    """ Greps output structure. """
    if not self.optimize: return self.input_structure

    from re import compile
    from itertools import chain
    from numpy import dot, array
    from numpy.linalg import inv
    from quantities import e, angstrom
    from ..crystal import Structure, into_cell
    from ..error import GrepError
    regex0 = compile('Final asymmetric unit coordinates') 
    regex1 = compile( 'Final fractional(?:/Cartesian)? coordinates '           \
                      'of (?:atoms|surface)' )
    with self.__stdout__() as file:
      # find beginning of final structure
      found = None
      for line in file:
        found = regex0.search(line)
        if found is not None: break
        found = regex1.search(line)
        if found is not None: break
      if found is None: raise GrepError('Could not grep atoms.')
      # find first atom
      regex = compile('(\d+)\s+([A-Z][a-z]?(?:\d+)?)\s+(?:c\s+)?\S(\s+\S+){3,5}')
      for line in file:
        found = regex.search(line.replace('*', ''))
        if found: break
      if found is None: raise GrepError('Could not grep atoms.')
      # create result
      result = Structure()
      result.cell = self.cell
      invcell = inv(self.cell)
      for line in chain([line], file):
        line = line.replace('*', '')
        found = regex.search(line)
        if found is None: break
        data = line.split()
        
        # figure out cartesian coords
        pos = array(data[3:6], dtype='float64')
        if self.dimensionality > 0: 
          d = self.dimensionality
          pos[:d] = dot(result.cell[:d,:d], pos[:d])
          pos[:d] = into_cell(pos, result.cell, invcell)[:d]

        # adds atom
        result.add_atom(pos=pos, type=data[1], label=int(data[0]))
        radius = float(data[6]) if len(data) > 6 else 0
        if abs(radius) > 1e-8: result[-1].radius = radius * angstrom
        charge = float(data[7]) if len(data) > 7 else 0
        if abs(charge) > 1e-8: result[-1].charge = charge * e
      # add more info
      try: result.energy = self.energy
      except: pass
      return result

  @property 
  @make_cached
  def input_structure(self):
    """ Greps input structure. """
    from re import compile
    from itertools import chain
    from numpy import dot, array
    from numpy.linalg import inv
    from quantities import e
    from ..crystal import Structure, into_cell
    from ..error import GrepError
    regex0 = compile('Fractional coordinates of asymmetric unit :')
    regex1 = compile('Mixed fractional/Cartesian coordinates of surface :')
    with self.__stdout__() as file:
      # find beginning of final structure
      found = None
      for line in file:
        found = regex0.search(line)
        if found is not None: break
        found = regex1.search(line)
        if found is not None: break
      if found is None: raise GrepError('Could not grep atoms.')
      # find first atom
      regex = compile('(\d+)\s+([A-Z][a-z]?(?:\d+)?)\s+\S(\s+\S+){5}')
      for line in file:
        found = regex.search(line)
        if found: break
      if found is None: raise GrepError('Could not grep atoms.')
      # create result
      result = Structure()
      result.cell = self.input_cell
      invcell = inv(result.cell)
      for line in chain([line], file):
        found = regex.search(line)
        if found is None: break
        data = line.replace('*', '').split()
        
        # figure out cartesian coords
        pos = array(data[3:6], dtype='float64')
        if self.dimensionality > 0: 
          d = self.dimensionality
          pos[:d] = dot(result.cell[:d,:d], pos[:d])
          pos[:d] = into_cell(pos, result.cell, invcell)[:d]

        # adds atom
        result.add_atom(pos=pos, type=data[1], label=int(data[0]))
        charge = float(data[6])
        if abs(charge) > 1e-8: result[-1].charge = charge * e
        occupancy = float(data[7])
        if abs(occupancy) > 1e-8: result[-1].occupancy = occupancy
      return result

  @property
  @make_cached
  def energy(self):
    """ Greps energy from result. """
    from re import M
    from quantities import eV
    from ..error import GrepError
    regex = 'Total\s+lattice\s+energy\s*=\s*(\S*)\s*eV'
    result = self._find_last_STDOUT(regex)
    if result is None:
      regex = 'Total\s+lattice\s+energy\s*:\s*\\n'                               \
              '\s*Primitive\s*unit\s*cell\s*=\s*(\S+)\s*eV'
      result = self._find_last_STDOUT(regex, M)
      if result is None: raise GrepError('Could not find energy.')
    return float(result.group(1)) * eV

  @property
  def _finished(self):
    """ True if job finished. """
    return self._find_first_STDOUT('Job Finished at') is not None

  @property
  def qeq_charges(self):
    """ Greps QEq charges from output. """
    from re import compile
    from numpy import array
    from quantities import e
    from ..error import GrepError
    if not self.qeq: raise GrepError('Not a QEq calculation.')
    with self.__stdout__() as file:
      found = False
      pattern = ['Final', 'charges', 'from', 'QEq', ':']
      for line in file:
        if line.split() == pattern: found = True; break
      if not found: raise GrepError('Could not find QEq charges.')

      pattern = compile('^\s*\d+\s+\d+\s+(\S*)\s*$')
      for line in file:
        if pattern.match(line) is not None: break
      result = [float(line.split()[-1])]
      for line in file:
        found = pattern.match(line)
        if found is None: break
        result.append(float(found.group(1)))
    return array(result) * e

class Extract(AbstractExtractBase, OutputSearchMixin, ExtractBase):
  """ Extracts GULP data from an output file. """
  def __init__(self, directory=None, **kwargs):
    """ Initializes extraction object. 
    
        :param directory: 
          Directory where the OUTCAR resides. 
          It may also be the path to an OUTCAR itself, if the file is not
          actually called OUTCAR.
    """
    from os.path import exists, isdir, basename, dirname
    from ..misc import RelativePath
       
    self.STDOUT = 'gulp.out'
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
    """ True if calculations ran to completion. """
    try: return self._finished
    except: return False
