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

class Eval:

  def __init__(self, size=10):
    from random import randint
    import numpy

    self.size = size
    self.target = numpy.array([randint(0,1) for u in xrange(self.size) ])
    self.nbevals = 0

  def __call__(self, indiv,comm = None):
    from math import fabs

    self.nbevals += 1
    result = 0e0
    for r in indiv.genes - self.target:
      result += fabs(r)
    return float(result) / float(len(indiv.genes))



    
def  main():
  import copy
  from boost.mpi import world
  import numpy
  from pylada.ga import darwin as dd, bitstring, standard, ce

  class Darwin: pass


  def stop_at_zero(self):
    for indiv in self.population:
      if indiv.fitness == 0e0: return False
    return True

  bitstring.Individual.size = 80
  darwin = Darwin()
  darwin.comm = world
  darwin.comm.do_print = darwin.comm.rank == 0

  evaluation = Eval()
  evaluation.target = numpy.array([1 for u in xrange(bitstring.Individual.size)])
  if darwin.comm.do_print: print "Target: ", evaluation.target

  darwin.evaluation = standard.mpi_population_evaluation( darwin, evaluation )

  def print_nb_eval(darwin):
    if not darwin.comm.do_print: return True
    print "Number of functional evaluations: ", evaluation.nbevals
    return True

  darwin.checkpoints = [ standard.print_offspring, 
                         standard.average_fitness,
                         standard.best, 
                         print_nb_eval,
                         stop_at_zero ]

  mating = standard.Mating(sequential=False)
  mating.add( bitstring.Crossover(rate=0.25), rate=0.8 )
  mating.add( bitstring.Mutation(rate=3e0/float(bitstring.Individual.size)), rate=0.2 )

  darwin.mating = standard.bound_method(darwin, standard.Mating(sequential=True))
  darwin.mating.add( mating, rate=0.8 )

  darwin.taboo = standard.bound_method(darwin, standard.Taboo(diversity=True))

  darwin.rate   = 0.1
  darwin.popsize = 100
  darwin.max_gen = 300

  dd.run(darwin)

if __name__ == "__main__":
  main()
