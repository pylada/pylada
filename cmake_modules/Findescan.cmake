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

if(ESCAN_INCLUDE_DIRS AND ESCAN_LIBRARY AND GENPOT_LIBRARY)
  set( ESCAN_FIND_QUIETLY TRUE)
endif(ESCAN_INCLUDE_DIRS AND ESCAN_LIBRARY AND GENPOT_LIBRARY)

find_path(_ESCAN_INCLUDE_DIRS
  NAMES
  coulomb_4pair_api.mod
  PATHS
  $ENV{ESCAN_INCLUDE_DIRS}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  escan nanopse
)
if(NOT _ESCAN_INCLUDE_DIRS)
  find_path(_ESCAN_INCLUDE_DIRS
    NAMES
    COULOMB_4PAIR_API.mod
    PATHS
    $ENV{ESCAN_INCLUDE_DIRS}
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
    escan nanopse
  )
endif(NOT _ESCAN_INCLUDE_DIRS)


FIND_LIBRARY(_ESCAN_LIBRARY
  pescan
  PATH
  $ENV{ESCAN_LIBRARY_DIR}
)
FIND_LIBRARY(_GENPOT_LIBRARY
  genpot
  PATH
  $ENV{ESCAN_LIBRARY_DIR}
)



INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(escan DEFAULT_MSG _ESCAN_LIBRARY _GENPOT_LIBRARY _ESCAN_INCLUDE_DIRS)
IF(_ESCAN_INCLUDE_DIRS AND _ESCAN_LIBRARY AND _GENPOT_LIBRARY)
  set(ESCAN_INCLUDE_DIRS ${_ESCAN_INCLUDE_DIRS} CACHE PATH "Path to nanopse include directory.")
  set(ESCAN_LIBRARY ${_ESCAN_LIBRARY} CACHE PATH "Path to escan library.")
  set(GENPOT_LIBRARY ${_GENPOT_LIBRARY} CACHE PATH "Path to genpot library.")
  SET(ESCAN_FOUND TRUE)
  unset(_ESCAN_INCLUDE_DIRS CACHE)
  unset(_GENPOT_LIBRARY CACHE)
  unset(_ESCAN_LIBRARY CACHE)
ENDIF(_ESCAN_INCLUDE_DIRS AND _ESCAN_LIBRARY AND _GENPOT_LIBRARY)
