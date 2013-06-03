set (ENUM_PY __init__.py transforms.py defects.py)
set (ENUM_CC module.cc ndimiterator.cc fciterator.cc manipulations.cc)
set (ENUM_H ndimiterator.h fciterator.h manipulations.h)

install(FILES ${ENUM_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/enum)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (pyenum SHARED ${ENUM_CC} ${ENUM_H})
set_target_properties(pyenum PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pyenum PROPERTIES PREFIX "")
set_target_properties(pyenum PROPERTIES SUFFIX ".so")
target_link_libraries(pyenum ${PYTHON_LIBRARIES})
install(TARGETS pyenum DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/enum)

add_subdirectory(tests)
