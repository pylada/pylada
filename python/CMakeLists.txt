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
