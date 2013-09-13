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

""" Subpackage containing extraction methods for CRYSTAL's property output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']

from ...tools import make_cached
from ...tools.extract import search_factory, AbstractExtractBase
from ...error import GrepError
OutputSearchMixin = search_factory('OutputSearchMixin', 'stdout', __name__, 'prop.out')

class ExtractBase(object):
  """ Implementation class for extracting CRYSTAL's property output. """
  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()

  @property
  @make_cached
  def is_bandstructure(self):
    """ True if contains band-structure calculation. """
    result = self._find_first_STDOUT("""^\s*\*\s*BAND\s+STRUCTURE\s*\*\s*$""")
    return result is not None

  @property
  @make_cached
  def end_date(self):
    """ Title of the calculations. """
    from datetime import datetime
    pattern = "E+\s+TERMINATION\s+DATE\s*(\d+)\s+(\d+)\s+(\d+)\s+TIME\s+(\S+)"
    regex = self._find_last_STDOUT(pattern)
    if regex is None: 
      raise GrepError( 'Could not grep end time from '                         \
                       '{0.directory}/{0.STDOUT}.'.format(self) )

    day    = int(regex.group(1))
    month  = int(regex.group(2))
    year   = int(regex.group(3))
    hour, minute, second = regex.group(4).split(':')
    hour   = int(hour)
    minute = int(minute)
    second = int(second[:second.find('.')])
    return datetime( year=year, month=month, day=day,
                     hour=hour, minute=minute, second=second )

  @property
  @make_cached
  def start_date(self):
    """ Title of the calculations. """
    from datetime import datetime
    search = " EEEEEEEEEE STARTING  DATE"
    with self.__stdout__() as file: 
      for line in file: 
        if len(line) <= len(search): continue
        if line[:len(search)] == search: break
      line = line.split()[3:]
      day, month, year = [int(u) for u in line[:3]]
      time = line[4].split(':')
      hour, minute, second = int(time[0]), int(time[1]), time[2]
      second = int(second[:second.find('.')])
      return datetime( year=year, month=month, day=day,
                       hour=hour, minute=minute, second=second )
  @property
  def spin_polarized(self):
    """ True if a spin-polarized calculation. """
    regex = self._find_first_STDOUT('SPIN POLARIZED')
    return regex is not None

  @property
  @make_cached
  def bandstructure(self):
    """ Computed band-structure. 

        If the band-structure was *not* computed, returns None.
        If it was, this is a named tuple ``(kpoints=[...], eigenvalues=[...])``.
        The kpoints are in the same order as per the actual calculation. They
        are always in fractional coordinates. None of that fake floating point
        iss thingie from crystal.  The eigenvalues are given in hartree, as a
        ``kpoint x band`` matrix.
    """
    from collections import namedtuple
    from re import M
    from os.path import join, exists
    from numpy import array
    from quantities import hartree
    from ... import CRYSTAL_propnames as propnames
    # if not a band-structure, then returns.
    if not self.is_bandstructure: return None
    bandfile = join(self.directory, propnames['BAND.DAT'].format('prop'))
    if not exists(bandfile):
      raise GrepError( 'Could not find bands-structure file {0}.'              \
                       .format(bandfile) )

    # create the result class.
    BandStructure = namedtuple('BandStructure', ['kpoints', 'eigenvalues'])

    # then greps k-points
    regex = r"^\s*LINE\s+\d+\s*\(\s*({0})\s+({0})\s+({0})\s*:\s*"              \
            r"\s*({0})\s+({0})\s+({0})\s*\)\s*({0})\s*POINTS\s*\-\s*"          \
            r"SHRINKING_FACTOR\s+({0})\s*$"                                    \
            r"\s*CARTESIAN\s+COORD\.\s+\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)"        \
            r"\s*:\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)\s*STEP\s+(\S+)\s*$"       \
            .format("(?:-|\+)?\d+")
    kpoints = []
    for kpre in self._search_STDOUT(regex, M):
      # if the segments are continuous, then CRYSTAL does not recompute the
      # first point of the next segment (since it is equal to the last point of
      # the current). We can figure that out by checking whether the computed
      # shrinking factor is a multiple of the number of points.
      # Finally, the same ratio tells us what the original shrinking factor is,
      # which helps us determine what the step is between two points.
      n, iss = int(kpre.group(7)), int(kpre.group(8))
      with_first = iss % n  != 0
      nbpoints = float(n-1 if with_first else n)
      shrink = float(iss // (n-1 if with_first else n))
      a = array([kpre.group(i) for i in xrange(1, 4)], dtype='float64')
      b = array([kpre.group(i) for i in xrange(4, 7)], dtype='float64')
      step = (b - a) / float(nbpoints) / float(shrink)
      first = a / shrink if with_first else a / shrink + step
      kpoints.extend( first + step * float(i) for i in xrange(n) )
    kpoints = array(kpoints)
    # If spin polarized calculation, then CRYSTAL print's kpoint stuff twice.
    # We don't want that.
    if self.spin_polarized: kpoints = kpoints.reshape(2, -1, 3)[0]

    # then greps eigenvalues from band-structure file.
    eigenvalues = []
    with open(bandfile, 'r') as file:
      for line in file: 
        line = line.rstrip()
        if line[0] != '@' and line[0] != '#':
          eigenvalues.append(line.split()[1:])

      eigenvalues = array(eigenvalues, dtype='float64') * hartree
    # If spin-polarized, then eigenvalues are spin x kpoints x band
    if self.spin_polarized:
      eigenvalues = eigenvalues.reshape(2, -1, eigenvalues.shape[-1])
    return BandStructure(kpoints, eigenvalues)

  @property
  @make_cached
  def structure(self):
    """ Structure as read by CRYSTAL. """
    from numpy import array
    from re import compile
    from ...crystal import Structure
    from ... import periodic_table as pt
    result = Structure()
    with self.__stdout__() as file:
      for line in file: 
        if "DIRECT LATTICE VECTOR COMPONENT" in line: break
      try: 
        result.cell = array( [file.next().split() for i in xrange(3)],
                             dtype='float64' ).T
      except StopIteration: 
        raise GrepError('Could not determine lattice parameter.')
      natoms = None
      for line in file: 
        if 'N. OF ATOMS PER CELL' in line:
          natoms = int(line.split()[5])
          break
      if natoms is None:
        raise GrepError('Could not determine number of atoms in unit-cell.')
      regex = compile( r"^\s*(\d+)\s+(\d+)\s+([A-Z]+)\s+\d+\s+"                \
                       r"(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s*$" )
      for line in file:
        found = regex.match(line)
        if found is None: continue
        Z = int(found.group(2))
        type = pt.find(atomic_number=Z%100).symbol
        result.add_atom( pos=[float(found.group(i)) for i in xrange(4, 7)],  
                         group=found.group(3),
                         index=int(found.group(1)),
                         Z=Z, type=type )
    return result


class Extract(AbstractExtractBase, OutputSearchMixin, ExtractBase):
  """ Extracts DFT data from an OUTCAR. """
  def __init__(self, directory=None, input=None, **kwargs):
    """ Initializes extraction object. 
    
        :param directory: 
          Directory where the OUTCAR resides. 
          It may also be the path to an OUTCAR itself, if the file is not
          actually called OUTCAR.
        :param input:
          Extraction object for self-consistent calculations. 
          If None, defaults to a normal crystal calculation in the input
          directory.
    """
    from os.path import exists, isdir, basename, dirname
    from ...misc import RelativePath
    from ..extract import Extract as ExtractSingle
       
    self.STDOUT = 'prop.out'
    """ Name of file to grep. """
    if directory is not None:
      directory = RelativePath(directory).path
      if exists(directory) and not isdir(directory):
        self.STDOUT = basename(directory)
        directory = dirname(directory)
    AbstractExtractBase.__init__(self, directory)
    ExtractBase.__init__(self)
    OutputSearchMixin.__init__(self)

    self.input = input if input is not None else ExtractSingle(self.directory)
    """ Extraction object for self-consistent calculations. """

  @property
  def success(self):
    from os.path import exists, join
    if not exists(join(self.directory, self.STDOUT)): return False
    try: self.end_date
    except: return False
    else: return True

