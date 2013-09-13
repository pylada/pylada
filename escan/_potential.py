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

""" Holds potential related stuff. """
__docformat__ = "restructuredtext en"
__all__ = [ 'AtomicPotential', 'localH', 'nonlocalH', 'soH' ]
from ..opt.decorators import broadcast_result

localH = 0
""" Local hamiltonian. """
nonlocalH = 1
""" Non-local hamiltonian. """
soH = 2
""" Spin-orbit hamiltonian. """

class AtomicPotential(object):
  """ Holds parameters to atomic potentials. """
  def __init__(self, path, nonlocal=None, s=None, p=None, d=None, pnl=None, dnl=None):
    from ..opt import RelativeDirectory

    self._filepath = RelativeDirectory(path=path)
    """ Private path to pseudopotential file. 
    
        Path is a relative directory for added transferability from computer to
        computer.
    """
    self._nonlocal = None if nonlocal is None else RelativeDirectory(nonlocal)
    """ Private path to non-local part, or None. 
    
        Path is a relative directory for added transferability from computer to
        computer.
    """
    self.s =  s if s is not None else 0
    """ s parameter """
    self.p =  p if p is not None else 0
    """ p parameter """
    self.d =  d if d is not None else 0
    """ d parameter """
    self.pnl = pnl if pnl is not None else 0
    """ p non-local parameter """
    self.dnl = dnl if dnl is not None else 0
    """ d non-local parameter """

  @property
  def filepath(self):
    """ Path to pseudopotential file. """
    return self._filepath.path
  @filepath.setter
  def filepath(self, value): self._filepath.path = value

  @property
  def nonlocal(self):
    """ Path to pseudopotential file. """
    return self._nonlocal.path
  @nonlocal.setter
  def nonlocal(self, value): self._nonlocal.path = value

  def __repr__(self):
    """ Tuple representation of self. """
    result = '"{0}"'.format(self._filepath.unexpanded)
    if self.nonlocal is None: result += ", None"
    else: result += ', "%s"' % (self._nonlocal.unexpanded)
    result += ", %f, %f, %f, %f, %f" % (self.s, self.p, self.d, self.pnl, self.dnl)
    return result

  @broadcast_result(key=True)
  def get_izz(self, comm = None):
    """ Returns izz string greped from pseudopotential file. """
    with open(self.filepath, "r") as file:
      return file.readline().split()[1]

