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

""" A GA subpackage defining standard genetic operator for bitstrings. """
__docformat__ = "restructuredtext en"

class VariableSizeCrossover(object):
  """ A crossover operation. """

  def __init__(self, nmin = -1, nmax = -1, step=2):
    """ Initializes crossover. 

        :Parameter:
          nmin : int
            Minimum size of  a bitstring. Only used if two bitstrings differ in
            size.
          nmin : int
            Minimum size of  a bitstring. Only used if two bitstrings differ in
            size.
    """
    self.nmin = nmin
    """ Minimum bitstring size. """
    self.nmax = nmax
    """ Maximum bitstring size. """

  def __call__(self, first, second):
    """ Create new individual from ``first`` and ``second`` using bitstring crossover. 
    
        :Parameter:
          first 
            First parent individual.
          second 
            This is the other parent individual.

        :return: An offspring individual.
    """
    from copy import deepcopy
    from random import randint
    from numpy import array

    a = [u for u in first.genes]
    b = [u for u in second.genes]
    if len(a) <= 1 : return deepcopy(second)
    if len(b) <= 1 : return deepcopy(first) 
    if len(a) <= 3 and len(b) <= 3: a = a[:1] + b[1:]
    elif len(a) <= 3: a += b[3:]
    elif len(b) <= 3: a = b + a[3:]
    else:
      i = randint(1, min(len(a), len(b))-2)
      j = randint(i+1, min(len(a), len(b)))
      a = a[:i] + b[i:j] + a[j:]
    result = deepcopy(first)
    result.genes  = array(a)
    assert len(result.genes) >= self.nmin 
    assert len(result.genes) <= self.nmax
    if hasattr(result, "fitness"): delattr(result, "fitness")

    return result
  
  def __repr__(self):
    return "{0.__class__.__name__}({0.nmin}, {0.nmax})".format(self)

class Crossover(object):
  """ A crossover operation. """

  def __init__(self, rate = 0.5):
    self.rate = rate

  def __call__(self, a, b):
    from random import uniform
    
    assert len(a.genes) == len(b.genes), "Bitstring must be of equal lengths"

    i = int( uniform(0, 1) * len(a.genes) )
    j = int( self.rate * len(a.genes) )
    if j >= len(a.genes): 
      a.genes[i:] = b.genes[i:] 
      a.genes[:j-len(a.genes)] = b.genes[:j-len(a.genes)]
    else: a.genes[i:j] = b.genes[i:j]

    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a

  def __repr__(self):
    return "{0.__class__.__name__}({0.rate})".format(self)

class Mutation(object):
  """ A mutation operation. """

  def __init__(self, rate = 0.10):
    """ Initializes mutation operation. 

        :Parameters:
          rate : float
            Mutation rate for each bit in the bitstring.
    """
    self.rate = rate
    """ Mutation rate for each bitstring. """

  def __call__(self, indiv):
    """ Returns invidual with mutated genes. """
    from random import uniform
    
    if len(indiv.genes) == 0: return indiv
    rate = self.rate if self.rate >= 0e0 else -float(self.rate)/float(len(indiv.genes))
    for i in xrange(len(indiv.genes)):
      if uniform(0, 1)  < rate:
        indiv.genes[i] = 0 if indiv.genes[i] == 1 else 1 
    if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
    return indiv

  def __repr__(self):
    return "{0.__class__.__name__}({0.rate})".format(self)
SwapMutation = Mutation
""" Fixed-size bitstring mutation operation. """


class GrowthMutation(object):
  """ Mutation which inserts a bit in the bitstring. """
  def __init__(self, nmin = -1, nmax = -1, step=1):
    """ Initiatizes the bit insertion.

        :Parameters:
          nmax : integer
            Maximum size of a bistring. If nmax = -1, there are no limits to the
            size of a bitstring.
          nmin : integer
            Minimum size of a bistring. If nmin = -1, there are no limits to the
            size of a bitstring.
    """
    self.nmax = nmax
    """ Maximum bistring size. """ 
    self.nmin = nmin
    """ Minimum bistring size. """ 
    assert self.nmin < self.nmax or self.nmin == -1 or self.nmax == -1,\
           ValueError("nmin and nmax are incorrect.")
    self.step = step
    """ By how much to grow or shrink. """

  def __call__(self, indiv):
    """ Inserts extra bit in bitstring. """
    from numpy import array
    from random import randint

    step = getattr(self, 'step', 1)
    l = [u for u in indiv.genes]
    if randint(0,1) == 0 and\
       ( (self.nmin == -1 and len(l) > step - 1) or \
         (self.nmin != -1 and len(l) > step - 1 + self.nmin) ):
      nmin = self.nmin if self.nmin != -1 else 0
      i = randint(nmin, len(l)-step)
      for j in xrange(step): l.pop(i)
    elif self.nmax == -1 or (len(l) <= self.nmax - step): 
      i = randint(0, len(l))
      j = randint(0, 1)
      for i in xrange(self.step): l.insert(i, j)
    indiv.genes = array(l)
    assert len(indiv.genes) >= self.nmin or self.nmin == -1
    assert len(indiv.genes) <= self.nmax or self.nmax == -1
     
    if hasattr(indiv, "fitness"): delattr(indiv, "fitness")
    return indiv

  def __repr__(self):
    return "{0.__class__.__name__}({0.nmin}, {0.nmax})".format(self)
