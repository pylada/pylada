if(do_tests)

  add_test(enum_cpp.py python ${CMAKE_CURRENT_SOURCE_DIR}/cpp.py)
  add_test(enum_transform.py python ${CMAKE_CURRENT_SOURCE_DIR}/transform.py)
  add_test(enum_supercell.py python ${CMAKE_CURRENT_SOURCE_DIR}/supercells.py)
  add_test(enum_generator.py python ${CMAKE_CURRENT_SOURCE_DIR}/generator.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(do_tests)
