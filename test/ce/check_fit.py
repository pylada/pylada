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

import numpy
from math import sqrt
import random
from pylada import crystal, ce
from pylada.mpi import world


lattice = crystal.bravais.fcc(); lattice.sites[0].type = "A", "B"
lattice.set_as_crystal_lattice()
mlclasses = ce.MLClusterClasses("input.xml", False)
fit = ce.PairRegulatedFit(mlclasses, lattice, alpha=3.0, tcoef=1, type="laks")

fit.read_directory("data")
A, b = fit()
x, residual, rank, s = numpy.linalg.lstsq(A, b)

ncls, nstr = fit.size()
size = max( int(float(nstr)*0.33333333333), 1 )
full_set = []
predsets = []
for i in xrange(80):
  if len(full_set) < size:
    full_set = [ j for j in range(nstr) ]
    random.shuffle(full_set)
  predsets.append( full_set[:size] )
  full_set = full_set[size:]
predsets = [
  [  8,  16,  15,  14,  17,  18,  10,  25 ], 
  [  1,  22,  24,   3,  21,  19,   5,   4 ], 
  [  6,   2,  12,  23,   9,  13,  20,   7 ], 
  [ 24,  10,   5,  11,  14,   2,  20,  18 ], 
  [  7,  13,   1,   4,  25,  12,  23,  22 ], 
  [ 19,  17,   8,   6,  15,  16,   9,  21 ], 
  [  2,   6,   1,   5,   8,  17,   9,  23 ], 
  [ 13,  10,   7,  25,  16,  18,   3,   4 ]] 
for i in xrange(len(predsets)):
  for j in xrange(len(predsets[i])):
    print "%3i " %(predsets[i][j] + 1),
    predsets[i][j] -= 1

errors = ce.leave_many_out(fit, predsets)
prediction = 0e0
training = 0e0
n = 0
for i in xrange(errors.shape[0]):
  npred, pred = 0, 0e0
  for j in xrange(errors.shape[1]):
    if j in predsets[i]:
      prediction += errors[i,j]*errors[i,j]
      pred += errors[i,j]*errors[i,j]
      n += 1
      npred += 1
    else: training += errors[i,j]*errors[i,j]
  print "set ", i, ": ", sqrt(pred / float(npred))
prediction /= float(n)
training /= float(errors.shape[0]*errors.shape[1] - n)
print errors.shape, n
print "Training error: ", sqrt(training)
print "Prediction  error: ", sqrt(prediction)
