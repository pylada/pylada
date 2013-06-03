set (PYCE_PY cluster.py __init__.py factory.py _single_site_factory.py)
set (CE_CC module.cc productilj.cc)
set (CE_H productilj.h)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ce)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyce SHARED ${CE_CC} ${CE_H})
set_target_properties(pyce PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyce PROPERTIES PREFIX "")
set_target_properties(pyce PROPERTIES SUFFIX ".so")
target_link_libraries(pyce ${PYTHON_LIBRARIES})
install(TARGETS pyce DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ce)


add_subdirectory(tests)
