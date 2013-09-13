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

""" A GA subpackage defining standard genetic operator for elemental alloys. """
__docformat__ = "restructuredtext en"
from operators import *
from converter import *
from functional import Darwin
from ...bitstring import Individual as BitstringIndividual

class Individual(BitstringIndividual):
  """ An individual for elemental superlattices. 
      
      Comparison between two individuals expect that a bitstring represents an
      elemental superlattice, with translation and inversion symmetries.
  """
  def __init__(self, nmin=0, nmax=20, step=2, dosym=False): 
    """ Initializes a bitstring individual randomly. """
    from random import randint
    import numpy

    assert nmin % 2 == 0 and nmax % 2 == 0 
    self.size = randint(nmin//2, nmax//2) * 2
    super(Individual, self).__init__(self.size)
    self.genes = numpy.array([ int(randint(0,1)) for i in xrange(self.size) ], dtype="int")
    self.dosym = dosym

  def __eq__(self, a): 
    """ Compares two elemental superlattices. """

    if a is None: return False
    if not hasattr(a, "genes"): return False
    if len(a.genes) != len(self.genes): return False
    if getattr(self, 'dosym', False): return all(a.genes - self.genes == 0) 
    dosym = getattr(self, 'dosym', 1)
    
    N = len(a.genes)
    for i in range(N):
      # periodicity
      found = True
      for j in range(N):
        if j % dosym != 0: continue
        if a.genes[ (i+j) % N ] != self.genes[j]:
          found = False
          break
      if found: return True

      # inversion + periodicity 
      if dosym != 1: continue
      found = True
      for j in range(N):
        if j % dosym != 0: continue
        if a.genes[ -((i+j) % N) ] != self.genes[j]:
          found = False
          break
      if found: return True

    return False
