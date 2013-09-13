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

if(eigen_INCLUDE_DIR)
  set( eigen_FIND_QUIETLY True)
endif(eigen_INCLUDE_DIR)

find_path(eigen_INCLUDE_DIR
  NAMES
  Eigen/Core
  Eigen/LU
  Eigen/Geometry
  Eigen/Cholesky
  PATHS
  $ENV{eigen_INCLUDE_DIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  eigen3 eigen2 eigen
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(eigen DEFAULT_MSG eigen_INCLUDE_DIR)

# gets wether eigen2 or eigen3
if(EIGEN_FOUND) 
  file(WRITE ${PROJECT_BINARY_DIR}/pylada_dummy.cc
        "#include <Eigen/Core>\n"  
        "#if not EIGEN_VERSION_AT_LEAST(3,0,0)\n"
        "#  error\n" 
        "#endif\n" )
  try_compile(is_eigen3  ${PROJECT_BINARY_DIR} ${PROJECT_BINARY_DIR}/pylada_dummy.cc
              COMPILE_DEFINITIONS -I${eigen_INCLUDE_DIR} 
              CMAKE_FLAGS -DCMAKE_CXX_LINK_EXECUTABLE="echo" 
              OUTPUT_VARIABLES _is_eigen3)
  file(REMOVE ${PROBJECT_BINARY_DIR}/pylada_dummy.cc)
endif(EIGEN_FOUND)
