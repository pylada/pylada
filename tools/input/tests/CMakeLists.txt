if(do_tests)

  add_test(functools_input_tree python ${CMAKE_CURRENT_SOURCE_DIR}/tree.py)
  add_test(functools_input_listblock python ${CMAKE_CURRENT_SOURCE_DIR}/listblock.py)
  add_test(functools_input_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)

endif(do_tests)
