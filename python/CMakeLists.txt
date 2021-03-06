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

set(PYLADA_PY onexit.py)

set(which_packages "'math', 'opt', 'crystal'") #, 'physics', 'periodic_table'")
if(pcm_functional) 
  set(which_packages "${which_packages}, 'pcm'")
endif(pcm_functional) 
if(ce_package) 
  set(which_packages "${which_packages}, 'ce'")
endif(ce_package) 
if(vasp_functional) 
  set(which_packages "${which_packages}, 'vasp'")
endif(vasp_functional) 
if(enum_package) 
  set(which_packages "${which_packages}, 'enum'")
endif(enum_package) 
if(jobs_package) 
  set(which_packages "${which_packages}, 'jobfolder'")
endif(jobs_package) 
if(ga_package) 
  set(which_packages "${which_packages}, 'ga'")
endif(ga_package) 
if(gulp_package) 
  set(which_packages "${which_packages}, 'gulp'")
endif(gulp_package) 
configure_file (
 "${PROJECT_SOURCE_DIR}/python/__init__.in.py"
 "${PROJECT_BINARY_DIR}/python/__init__.py"
)
install(FILES ${PROJECT_BINARY_DIR}/python/__init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)
install(FILES ${PYLADA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pypython SHARED module.cc)
set_target_properties(pypython PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pypython PROPERTIES PREFIX "")
set_target_properties(pypython PROPERTIES SUFFIX ".so")
target_link_libraries(pypython ${PYTHON_LIBRARIES})
install(TARGETS pypython DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)
install(FILES ${LIBCRYSTAL_PYSOURCES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)


set(PYHELPER_HEADERS python.h types.h wrap_numpy.h numpy_types.h
                     random_access_list_iterator.h
                     random_access_tuple_iterator.h
                     quantity.h
                     ppslot.hpp ppslot_crystal.hpp ppslot_math.hpp
                     ppslot_python.hpp )
set(PYHELPER_SOURCES object.cc quantity.cc)

if(NOT noinstall_headers)
  install(FILES ${PYHELPER_HEADERS} DESTINATION include/Pylada/python)
endif(NOT noinstall_headers)

add_subdirectory(tests)
