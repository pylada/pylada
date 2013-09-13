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

""" A GA subpackage to select figures for a cluster expansion. """
def _equivalent_sites(lattice):
  """ Returns a list containing only one site index per inequivalent sub-lattice. """
  from ...crystal import which_site

  result = set( i for i in range(len(lattice.sites)) ) 
  for i, site in enumerate(lattice.sites):
    if i not in result: continue
    for op in lattice.space_group:
      j = which_site( op(site.pos), lattice )
      if j != i and j in result: result.remove(j)
  return result

class Eval:

  def __init__( self, lattice, path="data", alpha=2, tcoef=10, lmo_ratio = 0.333333, **keywords):
    """ Creates an evaluation function for a multi-site cluster expansion.

        Requires on input the lattice, a path towards a directory of structure data,
        and the a short form for the allowed clusters.
    """
    from os.path import exists, join
    from re import compile, match
    from random import shuffle
    from ...ce import _create_clusters as Clusters, PairRegulatedFit

    self.lattice = lattice

    # data for the fit.
    self.alpha, self.tcoef = alpha, tcoef

    # checks that keywords are well formed.
    key_regex = compile("B(\d+)")
    for key in keywords.keys(): 
      a_re = match(key_regex, key)
      assert a_re is not None, "Keyword %s is not of the form B(\d+)" % (key)
      assert a_re.end() == len(key), "Keyword %s is not of the form B(\d+)" % (key)
      assert int(a_re.group(1)) > 1, "Cannot understand keyword %s" % (key)
      assert keywords[key] > 0, "Cannot understand input %s=%i" % (key, keywords[key])

    # sanity check.
    assert len(keywords) > 0, "No keywords found on input. Can't construct clusters.\n"

    keys = sorted(keywords.keys(), cmp)


    # now creates multi-lattice clusters, along with index bookkeeping.
    self._mlclasses = Clusters(lattice, nth_shell=0, order=0, site=0) # J0
    # creates J1 for these sites.
    for site in _equivalent_sites(self.lattice):
      self._mlclasses.extend(Clusters(lattice, nth_shell=0, order=1, site=site))
    # creates many bodies.
    self._fixed = len(self._mlclasses)
    for site in _equivalent_sites(self.lattice):
      for key in keys:
        regex = match(key_regex, key)
        order = int(regex.group(1))
        shell = keywords[key]
        self._mlclasses.extend(Clusters(lattice, nth_shell=shell, order=order, site=site))

    # creates fitting function 
    self.fitter = PairRegulatedFit(self._mlclasses, alpha=self.alpha, tcoef=self.tcoef)
    self.fitter.read_directory(path)

    # now reads/creates set.
    set_path = join(path, "set")
    self._sets = []
    ncls, nstr = self.fitter.size()
    if exists(set_path):
      # reads path 
      print "Found input set %s:\n" % (set_path)
      file = open(set_path, "r")
      for line in file:
        self._sets.append( [int(i) for i in line.split()] )
        print "  ", self._sets[-1]
        for i in self._sets[-1]:
          assert i < nstr, "Index in set is larger than the number of structures.\n"
      print
    else: 
      # creates set.
      size = max( int(float(nstr)*lmo_ratio), 1 )
      if size == 1: "Will perform leave-one out cross-validation.\n"
      else: 
        print "Creating leave-many out set:\n"
        full_set = []
        for i in xrange(size):
          if len(full_set) < size:
            full_set = [ j for j in range(nstr) ]
            shuffle(full_set)
          self._sets.append( full_set[:size] )
          full_set = full_set[size:]
          print "  ", self._sets[-1]

  def __len__(self):
    """ Returns the number of many-body figures (minus J0, J1, and pairs) """
    return len(self._mlclasses) - self._fixed
     


  def __call__(self, indiv):
    """ Returns cv score for this individual.
        Optimizes cv score with respect to pairs, and pair regularization.
    """
    from math import sqrt
    from numpy import arange, zeros, average
    from .. import ce

    assert len(self._mlclasses) - self._fixed == len(indiv.genes),\
           "Individual is of incorrect size.\n"
    onoffs = zeros((len(self._mlclasses),), dtype=bool)
    onoffs[:self._fixed] = True
    onoffs[self._fixed:] = indiv.genes

    self.fitter.extinguish_cluster(onoffs, mask=True) 
    self.fitter.alpha = self.alpha
    self.fitter.tcoef = self.tcoef

    if len(self._sets) == 0:
      errors = ce.leave_one_out(self.fitter)
      npreds, nstr = errors.shape
      prediction = errors[ arange(errors.shape[0]), arange(errors.shape[0]) ]
      return sqrt(average(prediction*prediction))
    else:
      errors = ce.leave_many_out(self.fitter, self._sets)
      npred, prediction = 0, 0e0
      for i in xrange(errors.shape[0]):
        for j in xrange(errors.shape[1]):
          if j not in self._sets[i]: continue
          prediction += errors[i,j] * errors[i,j]
          npred += 1
      return sqrt(prediction / float(npred))


class EvalFitPairs(Eval):
  """ Performs evaluation of many-body while optimizing regularised number of pairs. """

  max_nfuncs = 100
  max_pow = 15

  def __init__(self, pairs=(10, 20), *args, **kwarg):
    from math import fabs
    from ...ce import create_clusters as Clusters

    assert "B2" not in kwarg.keys(), "Cannot use B2 to initialize EvalFitPairs. Use pairs=(?,?).\n"

    if len(pairs) != 2 and len(pairs) !=3:
      raise ValueError, "Keyword pairs should be initialized with a 2 or 3-tuple.\n"
    if len(pairs) == 2: self.pairs = (pairs[0], pairs[1], 1)
    else:
      if fabs(pairs[0] -pairs[1]) < pairs[2]: 
        raise ValueError, "Third argument of keyword pairs is the step-size."
      self.pairs = (pairs[0], pairs[1], pairs[2]) 



    # creates self from base class.
    Eval.__init__(self, *args, **kwarg)
    # create pair terms.
    classes = None
    for site in _equivalent_sites(self.lattice):
      if classes is None:
        classes = Clusters(self.lattice, nth_shell=self.pairs[1], order=2, site=site)
        self._pairindex = [0]
        pairs = len(classes)
      else:
        self._pairindex.append(len(classes))
        classes.extend(Clusters(lself.attice, nth_shell=self.pairs[1], order=2, site=site))
        if len(classes) - self._pairindex[-1] > pairs: pairs = len(classes) - self._pairindex[-1]
    # adds pair terms.
    self.pairs = (self.pairs[0], pairs, self.pairs[2])
    self._pairindex.append(len(classes))
    self._fixed = len(classes) + 2 # pairs + J0 + J1 
    classes.extend(self._mlclasses)
    self._mlclasses = classes

    self.fitter.reset_clusterclasses(classes)


  def __call__(self, indiv):
    from numpy import arange, average, arrays, zeros
    from scipy.optimize import fmin as simplex
#   from scipy.optimize import fmin_powell as simplex
    from ...ce import leave_many_out, leave_one_out

    assert len(self._mlclasses) - self._fixed == len(indiv.genes),\
           "Individual is of incorrect size.\n"

    self.nfun = 0
    def callable_loo(x):
      from math import sqrt, log as ln, fabs, exp
      if self.nfun > EvalFitPairs.max_nfuncs: raise StopIteration
      if fabs(x[0]) > EvalFitPairs.max_pow: raise ValueError
      if fabs(x[1]) > EvalFitPairs.max_pow: raise ValueError
      self.nfun += 1
      self.fitter.alpha = x[0]
      self.fitter.tcoef = 100*exp(x[1])
      errors = leave_one_out(self.fitter)
      npreds, nstr = errors.shape
      prediction = errors[ arange(errors.shape[0]), arange(errors.shape[0]) ]
      return sqrt(average(prediction*prediction))
    def callable_lmo(x):
      from math import sqrt, log as ln, fabs, exp
      if self.nfun > EvalFitPairs.max_nfuncs: raise StopIteration
      if fabs(x[0]) > EvalFitPairs.max_pow: raise ValueError
      if fabs(x[1]) > EvalFitPairs.max_pow: raise ValueError
      self.nfun += 1
      self.fitter.alpha = x[0]
      self.fitter.tcoef = 100*exp(x[1])
      errors = leave_many_out(self.fitter, self._sets)
      npred, prediction = 0, 0e0
      for i in xrange(errors.shape[0]):
        for j in xrange(errors.shape[1]):
          if j not in self._sets[i]: continue
          prediction += errors[i,j] * errors[i,j]
          npred += 1
      return sqrt(prediction / float(npred))
    which_callable = callable_lmo
    if len(self._sets) == 0: which_callable = callable_loo

    # creates pair values.
    pairs = [ self.pairs[0] ]
    while pairs[-1] < self.pairs[1]:
      pairs.append( pairs[-1] + self.pairs[2] )
    if pairs[-1] != self.pairs[1]: pairs.append( self.pairs[1] )

    x0 = array([1.1,1.0])
    minvals = None
    for p in pairs:
      onoffs = zeros((len(self._mlclasses),), dtype=bool)
      for i, index in enumerate(self._pairindex[:-1]):
        nbpairs = min(i+p, self._pairindex[index+1])
        onoffs[i:nbpairs] = True
      onoffs[self._pairindex[-1]:self._fixed] = True
      onoffs[self._fixed:] = indiv.genes
      self.fitter.extinguish_cluster(onoffs, mask=True) 

#     x0, fopt, dummy, dummy, dummy, dummy = \
      fopt, iter, = 0,0
      try:
        x0, fopt, iter, dummy, dummy = \
            simplex(which_callable, x0, ftol=0.1, xtol=0.1, full_output=True, 
                    disp=False, maxfun=20, maxiter=20)
      except StopIteration:
        fopt = 1e12; x0 = array([1.1,1.0])
      except ValueError:
        fopt = 1e12; x0 = array([1.1,1.0])
      if minvals is None or minvals[0] > fopt:  minvals = (fopt, iter)

    return minvals[0]

class SetOver:
  """ Mating crossover operations over included sets of many-bodies.
      If inclusive, then the offspring contains all many-bodies that are in
      both parents, and some which are in either. If not inclusive, than the
      offspring contains some many-bodies which are set in either parents.
  """
  def __init__(self, rate=0.5, inclusive=True):
    """ Initializes SetOver parameters.
        inclusive species which type of \"SetOver\" to perform, and rate the
        acceptance rate of many-bodies. 
    """
    self.rate = rate; self.inclusive = inclusive

  def __call__(self, a, b): 
    from copy import deepcopy
    from random import random

    set_a = set()
    for i, n in enumerate(a.genes):
      if n: set_a.add(i)
    set_b = set()
    for i, n in enumerate(b.genes):
      if n: set_b.add(i)

    if self.inclusive: 
      a.genes[:] = False
      for i in set_a & set_b: 
        a.genes[i] = True
      for i in set_a ^ set_b: 
        if random() < self.rate:
          a.genes[i] = True
    else: 
      for i in set_a | set_b: 
        if random() < self.rate:
          a.genes[i] = True
    if hasattr(a, "fitness"): delattr(a, "fitness")
    return a;


class Taboo(object):
  """ Constrains individuals to having less than a given number of many-bodies. """
  def __init__(self, maxmbs = -1): self.maxmbs = maxmbs
  def __call__(self, darwin, indiv):
    """ Returns True if the individual has more than self.maxmbs manybodies. """
    if self.maxmbs <= 0: return False

    sum = 0
    for i in indiv.genes:
      if i: sum += 1 
    return sum > self.maxmbs 


    

class Individual(object):
  """ An individual for boolean bitstrings.

      The size of the individual is given statically by Individual.size
      @attention: This class does handle mpi at all. Created instances will
        differ from one process to the next.
  """

  size = 10
  max_mbs_oninit = -1

  def __init__(self, darwin):
    """ Initializes a boolean bitstring individual randomly.

        @attention: This class does handle mpi at all. Created instances will
          differ from one process to the next.
    """
    from random import choice, shuffle, randint
    from numpy import array

    super(Individual, self).__init__()
    self.genes = array([ choice([True,False]) for i in xrange(Individual.size) ], dtype=bool)

    if self.max_mbs_oninit > 0:
      list_ = [i for i in xrange(len(self.genes))]
      shuffle(list_)
      self.genes[:] = False
      self.genes[list_[:randint(1, min(20, Individual.size))]] = True
  
  def __eq__(self, a): 
    from math import fabs
    if a is None: return False
    for i, b in enumerate(a.genes):
      if self.genes[i] != b: return False
    return True
  
  def __str__(self):
    string = ""
    for i in self.genes:
      if i: string += "1"
      else: string += "0"
    return string
