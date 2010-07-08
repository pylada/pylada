if(GSL_INCLUDE_DIRS AND GSL_LIBRARY AND GSLCBLAS_LIBRARY)
  set( GSL_FIND_QUIETLY True)
endif(GSL_INCLUDE_DIRS AND GSL_LIBRARY AND GSLCBLAS_LIBRARY)

find_path(_GSL_INCLUDE_DIRS
  NAMES
  gsl_vector.h
  gsl_matrix.h
  PATHS
  $ENV{GSL_INCLUDE_DIRS}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  gsl
)

FIND_LIBRARY(_GSL_LIBRARY
  gsl
  PATH
  $ENV{GSL_LIBRARY_DIR}
)

FIND_LIBRARY(_GSLCBLAS_LIBRARY
  gslcblas
  PATH
  $ENV{GSL_LIBRARY_DIR}
)


IF(_GSL_INCLUDE_DIRS AND _GSL_LIBRARY AND _GSLCBLAS_LIBRARY)
  SET( GSL_INCLUDE_DIRS ${_GSL_INCLUDE_DIRS} CACHE
       PATH "Path to Gnu Scientific Librarie's includes." )
  unset(_GSL_INCLUDE_DIRS CACHE)
  SET( GSL_LIBRARY ${_GSL_LIBRARY} CACHE
       PATH "Path to Gnu Scientific Librarie's main library." )
  unset(_GSL_LIBRARY CACHE)
  SET( GSLCBLAS_LIBRARY ${_GSLCBLAS_LIBRARY} CACHE
       PATH "Path to Gnu Scientific Librarie's main library." )
  unset(_GSLCBLAS_LIBRARY CACHE)
  SET(GSL_FOUND TRUE CACHE BOOL "Whether the GSL libraries have been found.")
ENDIF(_GSL_INCLUDE_DIRS AND _GSL_LIBRARY AND _GSLCBLAS_LIBRARY)
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(gsl DEFAULT_MSG GSL_LIBRARY GSL_INCLUDE_DIRS)
