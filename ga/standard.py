#
#  Version: $Id$
#
""" Holds standard genetic algorithm operations. """
def no_taboo( self, _indiv ):
  """ Does nothing. """
  return False
def taboo( self, _indiv ):
  """ taboo makes sure that no two individuals in the population and the offspring are the same. """
  return _indiv in self.population or _indiv in self.offspring

def tournament( self, size = 2 ):
  """ deterministic tournament """
  import random
  list_ = range(len(self.population))
  random.shuffle(list_)
  list_ = list_[:size]
  result = 0
  for b in list_:
    if self.cmp_indiv(self.population[b], self.population[result]):
      result = b;
  return result

def cmp_indiv( a, b, tolerance = 1e-12 ):
  """ Compares two individuals. """
  from math import fabs
  if fabs(a.fitness - b.fitness) <  tolerance: return 0
  elif a.fitness < b.fitness: return -1
  return 1

def average_fitness(self):
  """ Prints out average fitness. """
  result = 0e0
  for indiv in self.population:
    result += indiv.fitness
  print "  Average Fitness: ", result / float(len(self.population))
  return True

def best(self):
  best = None
  for indiv in self.population:
    if best == None or  self.cmp_indiv( best, indiv ) == 1: 
      best = indiv
  print "  Best Individual: ", best, best.fitness
  return True

def print_population(self):
  print "  Population: "
  for indiv in self.population:
    print "    ", indiv, indiv.fitness
  return True

def print_offspring(self):
  print "  Offspring: "
  for indiv in self.population:
    if indiv.birth == self.current_gen - 1: 
      print "    ", indiv, indiv.fitness
  return True

def check_generation( self ):
  """ returns false if maximum number of generations was passed. """
  if self.max_gen < 0: return True
  return self.current_gen < self.max_gen
  
def add_population_evaluation(self, evaluation):
  """ Standard population evaluation. 
      Evaluates individual only if fitness attribute does not exist. 
      Fitness is the return of evaluation subroutine given on input.
      evaluation subroutine should take an individual at its only argument.
  """
  def popeval(self):
    for indiv in self.population:
      if not hasattr(indiv, "fitness" ): 
        indiv.fitness = self.indiv_evaluation(indiv)
    for indiv in self.offspring:
      if not hasattr(indiv, "fitness" ): 
        indiv.fitness = self.indiv_evaluation(indiv)
  self.indiv_evaluation = evaluation
  self.evaluation = popeval
  return self


class Mating(object):
  """ Aggregator of mating operators. 
      Mating operations can be added using the Mating.add bound method.
      Operations are called either sequentially[==and] or proportionally[==or] (either all
      or only one can be called). 
  """
      

  def __init__(self, sequential=False): 

    self.operators = []
    self.sequential = sequential

  def add(self, function, rate=1):


    # checks for number of arguments to function.
    nb_args = -1
    if hasattr(function, "__class__") and issubclass(function.__class__, Mating): pass
    elif hasattr(function, "func_code"): nb_args = function.func_code.co_argcount
    else: nb_args = function.__call__.im_func.func_code.co_argcount - 1

    if rate <= 0e0: raise ValueError, "rate argument cannot be negative (%s)." % (rate)
    self.operators.append( (function, rate, nb_args) )
  
  def __call__(self, darwin):

    from random import random

    def call_function(function, n, indiv = None):
      """ Calls functions/functors mating operations. """
      from copy import deepcopy
      # calls other Mating instances.
      if n == -1: return function(darwin)

      individuals = []
      if indiv != None: individuals.append(indiv)
      else: individuals.append( deepcopy(darwin.population[darwin.selection(darwin)]) )

      # calls unaries
      if   n == 1: return function( individuals[0] )

      # calls binaries
      b = individuals[0]
      while( b not in individuals ): b = darwin.population[darwin.selection(darwin)]
      individuals.append(b)
      if n == 2: return function( individuals[0], individuals[1] )

      # calls ternaries
      b = individuals[0]
      while( b not in individuals ): b = darwin.population[darwin.selection(darwin)]
      individuals.append(b)
      if n == 3: return function( individuals[0], individuals[1], individuals[2] )

      raise "Mating operations has to be unary, binary, or ternary.\n"

    indiv = None
    if self.sequential: # choose any operator depending on rate.
      while indiv == None: # makes sure we don't bypass all mating operations
        for function, rate, n in self.operators:
          if random() < rate: indiv = call_function( function, n, indiv )
    else: # choose only one operator.
      max = 0e0
      for function, rate, n in self.operators: max += rate
      assert rate > 0e0;
      last = 0e0
      r = random() * max
      for function, rate, n in self.operators:
        if r <= last + rate:
          indiv = call_function( function, n )
          break
        last += rate

    assert indiv != None, "%s" % (self.sequential)
    return indiv


def add_checkpoint(self, _chk):
  """ Adds a checkpoint """
  try: self.checkpoints.append( _chk ) 
  except AttributeError: self.checkpoints = [_chk]
  return self;

def fill_attributes(self):
  """ Checks self for correct attributes.
      Fills in where possible:
        _ "taboo" defaults to standard.taboo
        _ "selection" defaults to standard.selection
        _ "population" defaults to empty list []
        _ "popsize" defaults to len(self.population)
        _ "offspring" defaults to empty list []
        _ "rate" defaults to len(self.offspring)/self.popsize
        _ standard.check_generation is ALWAYS added to the checkpoints
        _ "current_gen" defaults to 0 
        _ "max_gen" defaults to 100, or current_gen+100 if max_gen > current_gen.
        _ "cmp_indiv" defaults to standard.cmp_indiv
  """
  import darwin 
  # must have an evaluation function.
  assert hasattr(self, "evaluation"), "No evaluation function!" 

  # Checks that self has an object Individual
  if not hasattr(self, "Individual"):
    from bitstring import Individual
    self.Individual = Individual

  # Checks whether self has a taboo object.
  if not hasattr(self, "taboo"): self.taboo = taboo

  # Checks whether self has a selection object.
  if not hasattr(self, "selection"): self.selection = tournament

  # checks whether there is a population.
  if not hasattr(self, "population"): self.population = []
  
  # checks whether there is a popsize.
  if not hasattr(self, "popsize"): self.popsize = len(self.population)

  # makes sure we have something to do.
  assert self.popsize != 0, "population or popsize attributes required on input."

  # checks whether there is a population.
  if not hasattr(self, "offspring"): self.offspring = []
  
  # checks whether there is a popsize.
  if not hasattr(self, "rate"):
    self.rate = float(len(self.offspring)) / float(self.popsize)

  # makes sure we have something to do.
  assert self.rate > float(0), "offspring or rate attributes required on input."
 
  # checks whether there is a checkpoint.
  self = add_checkpoint(self, check_generation)

  # checks current generation.
  if not hasattr(self, "current_gen"): self.current_gen = 0
  elif not hasattr(self, "max_gen"): self.max_gen = 100
  elif self.max_gen < self.current_gen: self.max_gen += 100

  if not hasattr(self, "max_gen"): self.max_gen = self.current_gen + 100

  if not hasattr(self, "cmp_indiv"): self.cmp_indiv = cmp_indiv
  return self
