include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pymath SHARED module.cc noopt.cc)
SET_SOURCE_FILES_PROPERTIES(noopt.cc PROPERTIES COMPILE_FLAGS -O0)
set_target_properties(pymath PROPERTIES OUTPUT_NAME math)
set_target_properties(pymath PROPERTIES PREFIX "")
set_target_properties(pymath PROPERTIES SUFFIX ".so")
target_link_libraries(pymath ${PYTHON_LIBRARIES})
install(TARGETS pymath DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)

set(PYMATH_HEADERS math.h eigen.h exceptions.h fuzzy.h misc.h
                   symmetry_operator.h algorithms.h )
set(PYMATH_SOURCES gruber.cc smith_normal_form.cc methods.cc)

if(NOT noinstall_headers)
  install(FILES ${PYMATH_HEADERS} DESTINATION include/Pylada/math)
endif(NOT noinstall_headers)

add_subdirectory(tests)
