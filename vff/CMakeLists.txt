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

set( PYVFF_HEADERS 
     node/pybase.h node/sequence.hpp node/iterator.hpp
     edge/pybase.h )
set( PYVFF_EXTRA_SOURCES  zb.cc )
set( PYVFF_SOURCES 
     module.cc
     node/pybase.cc 
     edge/pybase.cc )
set( PYVFF_PY __init__.py functional.py vff.py extract.py )

install(FILES ${PYVFF_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vff)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyvff SHARED ${PYVFF_SOURCES} ${PYVFF_HEADERS})
set_target_properties(pyvff PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyvff PROPERTIES PREFIX "")
set_target_properties(pyvff PROPERTIES SUFFIX ".so")
target_link_libraries(pyvff ${PYTHON_LIBRARIES})
install(TARGETS pyvff DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vff)


add_subdirectory(tests)
