if(do_tests)

  add_test(gulp_extract_structure python ${CMAKE_CURRENT_SOURCE_DIR}/extract_structure.py)
  add_test(gulp_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)
  add_test(gulp_functional python ${CMAKE_CURRENT_SOURCE_DIR}/functional.py)

  if(do_external_program_tests)

    add_test(gulp_run python ${CMAKE_CURRENT_SOURCE_DIR}/run.py)

  endif(do_external_program_tests)

endif(do_tests)
