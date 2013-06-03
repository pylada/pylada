CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (Pylada)
set (Pylada_VERSION_MAJOR 5)
set (Pylada_VERSION_MINOR 0)

option(noinstall_headers    "Will not install header files." on)
mark_as_advanced(noinstall_headers)

# Project package options 
option(ewald_functional       "Include Point Charge Ion functional."       off)
option(vasp_functional        "Compile vasp functional wrapper."           on)
option(crystal_functional     "Install CRYSTAL functional wrapper."        on)
option(vff_functional         "Compile Valence Force Field method."        on)
option(gulp_functional        "Include GULP wrappers."                     on)
option(ce_package             "Include Cluster Expansion package."         on)
option(enum_package           "Include supercell enumeration package"      on)
option(jobfolder_package      "Installs jobfolder python module."          on)
option(ipython_interface      "Include IPython interface."                 on)
option(ga_package             "Include genetic algorithm package"          on)
option(do_tests               "Compile tests."                             on)
option(do_external_program_tests "Allow tests which run external programs (vasp, crystal...)" on)

if(do_tests)
  ENABLE_TESTING()
endif(do_tests)

# Location where cmake first looks for cmake modules.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# Location where cmake looks for libraries
# adds environment SOFTWARE_ROOT if it exists.
# adds SOFTWARE_ROOT if it exists.
# adds environment LD_LIBRARY_PATH if it exists.
if (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
  set(CMAKE_LIBRARY_PATH $ENV{SOFTWARE_ROOT}/lib $ENV{SOFTWARE_ROOT}/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH $ENV{SOFTWARE_ROOT}/include)
else (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
  set(CMAKE_LIBRARY_PATH $ENV{HOME}/usr/lib $ENV{HOME}/usr/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH $ENV{HOME}/usr/include)
endif (NOT "$ENV{SOFTWARE_ROOT}" STREQUAL "")
if(SOFTWARE_ROOT)
  set(CMAKE_LIBRARY_PATH ${SOFTWARE_ROOT}/lib ${SOFTWARE_ROOT}/lib64  ${CMAKE_LIBRARY_PATH})
  set(CMAKE_INCLUDE_PATH ${SOFTWARE_ROOT}/include)
endif (SOFTWARE_ROOT)
if (NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
  string(REPLACE ":" ";" _LDLIBPATH $ENV{LD_LIBRARY_PATH})
  set(CMAKE_LIBRARY_PATH ${_LDLIBPATH} ${CMAKE_LIBRARY_PATH})
endif (NOT "$ENV{LD_LIBRARY_PATH}" STREQUAL "")
  


# Languages
enable_language(C)
enable_language(CXX)
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE "FC_" )
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)

# look for python version.
if(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)
  execute_process( 
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print \"%i.%i\" % sys.version_info[:2]"
    OUTPUT_VARIABLE PYTHON_VERSION
  )
  if( PYTHON_VERSION )
    string (STRIP ${PYTHON_VERSION} PYTHON_VERSION)
    set(PYTHON_VERSION ${PYTHON_VERSION} CACHE STRING "Version of the Python interpreter.")
    mark_as_advanced(PYTHON_VERSION)
    MESSAGE(STATUS "[Python] Version: ${PYTHON_VERSION}")
  else( PYTHON_VERSION )
    MESSAGE(STATUS "Could not determine python version.")
  endif( PYTHON_VERSION )
endif(NOT PYTHON_VERSION AND PYTHONINTERP_FOUND)
# sets python site-packages
if(NOT CMAKE_PYINSTALL_PREFIX)
  if(PYTHON_VERSION)
    set( CMAKE_PYINSTALL_PREFIX lib64/python${PYTHON_VERSION}/site-packages/ 
         CACHE PATH "Python modules will be installed here." )
  else(PYTHON_VERSION)
    execute_process( 
      COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
      OUTPUT_VARIABLE CMAKE_PYINSTALL_PREFIX
    )
    if( CMAKE_PYINSTALL_PREFIX )
      string (STRIP ${CMAKE_PYINSTALL_PREFIX} CMAKE_PYINSTALL_PREFIX)
      set(CMAKE_PYINSTALL_PREFIX ${CMAKE_PYINSTALL_PREFIX} CACHE STRING "Version of the Python interpreter.")
      mark_as_advanced(CMAKE_PYINSTALL_PREFIX)
      MESSAGE(STATUS "[Python] Version: ${CMAKE_PYINSTALL_PREFIX}")
    endif( CMAKE_PYINSTALL_PREFIX )
  endif(PYTHON_VERSION)
endif(NOT CMAKE_PYINSTALL_PREFIX)


# required packages
#   -- boost packages
if(ewald_functional)
  find_package(Boost 1.33.0 REQUIRED COMPONENTS math_c99 math_c99l math_c99f )
else(ewald_functional)
  find_package(Boost 1.33.0 REQUIRED)
endif(ewald_functional)
include_directories(${Boost_INCLUDE_DIRS})
#   -- eigen2
find_package( eigen REQUIRED )
include_directories(${eigen_INCLUDE_DIR})
set(PYLADA_WITH_EIGEN3 ${is_eigen3})

find_package(MPI) # only required for some tests in process and jobfolders
FortranCInterface_VERIFY()
find_package(numpy REQUIRED) # has to be after FortranCInterface_VERIFY for some reason

#include root source directory
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_BINARY_DIR}")

# the most basic stuff to compile/include
# add_subdirectory(opt)
add_subdirectory(python)
add_subdirectory(math)
add_subdirectory(crystal)
add_subdirectory(physics)
add_subdirectory(config)
add_subdirectory(tools)
add_subdirectory(misc)
add_subdirectory(process)
add_subdirectory(periodic_table)
add_subdirectory(errors)
if(ga_package)
  add_subdirectory(ga)
endif(ga_package)

if(ipython_interface)
  set(jobfolder_package on)
endif(ipython_interface)

if( vasp_functional )
  set(DOVASP True)
  add_subdirectory(vasp)
endif( vasp_functional )

if( crystal_functional )
  set(DOCRYSTAL True)
  add_subdirectory(dftcrystal)
endif( crystal_functional )

if( vff_functional )
  add_subdirectory(vff)
endif( vff_functional )

if( gulp_functional )
  add_subdirectory(gulp)
endif( gulp_functional )

if( ce_package )
  add_subdirectory(ce)
  set(enum_package True)
endif( ce_package )
if( enum_package )
  add_subdirectory(enum)
endif( enum_package )

if(ewald_functional)
  add_subdirectory(ewald)
endif(ewald_functional)

if(ipython_interface)
  add_subdirectory(ipython)
endif(ipython_interface)

if(jobfolder_package)
  add_subdirectory(jobfolder)
endif(jobfolder_package)

# Configure file.
if (CMAKE_BUILD_TYPE)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(PYLADA_DEBUG TRUE)
  endif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
endif (CMAKE_BUILD_TYPE)
configure_file (
 "${PROJECT_SOURCE_DIR}/PyladaConfig.h.in"
 "${PROJECT_BINARY_DIR}/PyladaConfig.h"
)

# adds root_exception to install.
if(NOT noinstall_headers)
  install(FILES root_exceptions.h DESTINATION include/Pylada)
endif(NOT noinstall_headers)
