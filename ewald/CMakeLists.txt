set(PYEWALD_PY __init__.py)
set(PYEWALD_F90 ep_com.f90 ewaldf.f90)
set(PYEWALD_HEADERS ewald.h)
set(PYEWALD_SOURCES ewald.cc module.cc erfc.cc)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyewald SHARED ${PYEWALD_SOURCES} ${PYEWALD_HEADERS} ${PYEWALD_F90})
set_target_properties(pyewald PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyewald PROPERTIES PREFIX "")
set_target_properties(pyewald PROPERTIES SUFFIX ".so")
target_link_libraries(pyewald ${PYTHON_LIBRARIES})
target_link_libraries(pyewald ${Boost_MATH_C99_LIBRARY}
                              ${Boost_MATH_C99L_LIBRARYl}
                              ${BOOST_C99F_LIBRARY} )
install(TARGETS pyewald DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ewald)
install(FILES ${PYEWALD_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ewald)

if(do_tests)
  add_subdirectory(tests)
endif(do_tests)
