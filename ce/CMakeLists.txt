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

set (PYCE_PY cluster.py __init__.py factory.py _single_site_factory.py)
set (CE_CC module.cc productilj.cc)
set (CE_H productilj.h)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ce)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyce SHARED ${CE_CC} ${CE_H})
set_target_properties(pyce PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyce PROPERTIES PREFIX "")
set_target_properties(pyce PROPERTIES SUFFIX ".so")
target_link_libraries(pyce ${PYTHON_LIBRARIES})
install(TARGETS pyce DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ce)


add_subdirectory(tests)
