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
