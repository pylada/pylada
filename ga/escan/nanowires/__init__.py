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

""" A GA subpackage defining standard genetic operator for nanowires. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Individual']

from ...bitstring import Individual as BitstringIndividual
class Individual(BitstringIndividual):
  """ An individual for nanowires.
      
      Comparison between two individuals expect that a bitstring represents an
      elemental nanowires.
  """
  def __init__(self, nmin=0, nmax=20):
    """ Initializes a bitstring individual randomly. """
    from random import randint
    import numpy

    self.size = randint(nmin, nmax)
    super(Individual, self).__init__(self.size)
    self.genes = numpy.array([ int(randint(0,1)) for i in xrange(self.size) ], dtype="int")
  
  def __eq__(self, a): 
    """ Compares two elemental nanowires. """
    from numpy import all
    if a is None: return False
    if not hasattr(a, "genes"): return False
    if len(a.genes) != len(self.genes): return False
    return all(self.genes == a.genes)

def exec_input(script, namespace = None):
  """ Executes an input script including namespace for ga nanowires. """ 
  from ....escan import exec_input as escan_exec_input
  from functional import Darwin
  from converter import Converter
  from ..evaluator import Bandgap as BandGapEvaluator, Dipole as DipoleEvaluator

  dictionary = {}
  if namespace is not None: dictionary.update(namespace)
  dictionary['Individual'] = Individual
  dictionary['Darwin'] = Darwin
  dictionary['Converter'] = Converter
  dictionary['BandGapEvaluator'] = BandGapEvaluator
  dictionary['DipoleEvaluator']  = DipoleEvaluator
  return escan_exec_input(script, dictionary)

def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for ga nanowires. """ 
  from ....escan import read_input as escan_read_input
  from functional import Darwin
  from converter import Converter
  from ..evaluator import Bandgap as BandGapEvaluator, Dipole as DipoleEvaluator

  dictionary = {}
  if namespace is not None: dictionary.update(namespace)
  dictionary['Individual'] = Individual
  dictionary['Darwin'] = Darwin
  dictionary['Converter'] = Converter
  dictionary['BandGapEvaluator'] = BandGapEvaluator
  dictionary['DipoleEvaluator']  = DipoleEvaluator
  return escan_read_input(filepath, dictionary)
