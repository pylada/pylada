if(do_tests)

  set(TYPES 0 1)
  foreach(TYPE ${TYPES})
    add_executable(fuzzy${TYPE} fuzzy.cc)
    set_property( TARGET fuzzy${TYPE} PROPERTY COMPILE_DEFINITIONS
                  PYLADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzy${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzy${TYPE})
  endforeach(TYPE)

  set(TYPES 0 1 2 3)
  foreach(TYPE ${TYPES})
    add_executable(fuzzyeigen${TYPE} fuzzyeigen.cc)
    set_property( TARGET fuzzyeigen${TYPE} PROPERTY COMPILE_DEFINITIONS
                  PYLADA_TEST_INCTYPE=${TYPE})
    ADD_TEST(math_fuzzyeigen${TYPE} ${CMAKE_CURRENT_BINARY_DIR}/fuzzyeigen${TYPE})
  endforeach(TYPE)

  add_executable(symtests symtests.cc)
  ADD_TEST( math_symtests ${CMAKE_CURRENT_BINARY_DIR}/symtests)

  add_library (gruber SHARED "gruber.cc")
  target_link_libraries(gruber ${PYTHON_LIBRARIES})
  set_target_properties(gruber PROPERTIES PREFIX "_")
  set_target_properties(gruber PROPERTIES SUFFIX ".so")
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS PYLADA_LIM=5)
  else("${CMAKE_BUILD_TYPE}" STREQUAL "Release") 
    set_property(TARGET gruber PROPERTY COMPILE_DEFINITIONS PYLADA_LIM=2)
  endif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

  add_library (smith SHARED "smith_form.cc")
  set_target_properties(smith PROPERTIES PREFIX "_")
  set_target_properties(smith PROPERTIES OUTPUT_NAME "smith")
  set_target_properties(smith PROPERTIES SUFFIX ".so")
  target_link_libraries(smith ${PYTHON_LIBRARIES})

  add_test(math_basics python ${CMAKE_CURRENT_SOURCE_DIR}/basics.py)
  add_test(math_pygruber python ${CMAKE_CURRENT_SOURCE_DIR}/gruber.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(math_pysmith python ${CMAKE_CURRENT_SOURCE_DIR}/smith.py ${CMAKE_CURRENT_BINARY_DIR})

endif(do_tests)
