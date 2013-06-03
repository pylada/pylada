set(LIBCRYSTAL_PYSOURCES __init__.py read.py write.py ABX.py A2BX4.py
                         bravais.py binary.py iterator.py)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})

add_library (pycrystal SHARED module.cc)
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
set_target_properties(pycrystal PROPERTIES SUFFIX ".so")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/crystal)
install(FILES ${LIBCRYSTAL_PYSOURCES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/crystal)
target_link_libraries(pycrystal ${PYTHON_LIBRARIES})

set(PYCRYSTAL_SOURCES atom/pybase.cc atom/cdi.cc atom/members.cc atom/getset.cc
                      structure/pybase.cc structure/cdi.cc structure/members.cc
                      structure/getset.cc structure/sequence.cc
                      structure/iterator.cc 
                      hart-forcade/pybase.cc hart-forcade/cdi.cc
                      hart-forcade/get.cc hart-forcade/members.cc
                      utilities.cc map_sites.cc equivalent_structures.cc
                      primitive.cc space_group.cc neighbors.cc
                      coordination_shells.cc confsplit.cc periodic_dnc.cc
                      methods.cc )


set(PYCRYSTAL_HEADERS atom/atom.h
                      structure/structure.h
                      hart-forcade/hart-forcade.h hart-forcade/macro.hpp
                      utilities.h supercell.h algorithms.h )

if(NOT noinstall_headers)
  install(FILES ${PYCRYSTAL_HEADERS} DESTINATION include/Pylada/crystal)
endif(NOT noinstall_headers)

add_subdirectory(tests)
