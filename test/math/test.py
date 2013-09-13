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

from numpy import array
from pylada.math import check_extract, check_frompy_to_cpp_vec, check_frompy_to_cpp_mat,\
                      check_topy_from_cpp_mat, check_topy_from_cpp_vec

mat = array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
vec = array([-1,0.5,0.5])

check_extract(mat)
print "extraction - OK"

check_frompy_to_cpp_mat(mat)
print "automatic matrix extraction - OK"
check_frompy_to_cpp_vec(vec)
print "automatic vector extraction - OK"

print check_topy_from_cpp_mat()
print "automatic matrix conversion - OK"
print check_topy_from_cpp_vec()
print "automatic vector conversion - OK"

