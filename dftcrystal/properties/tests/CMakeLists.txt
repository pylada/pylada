if(do_tests)
  # All will require external program to run tests.
  # This is because even the most basic operation on crystal structures runs the
  # program in the background
  if(do_external_program_tests)

    add_test(dftcrystal_properties_keywords python ${CMAKE_CURRENT_SOURCE_DIR}/keywords.py)
    add_test(dftcrystal_properties_bandstructure python ${CMAKE_CURRENT_SOURCE_DIR}/bandstructure.py)
    add_test(dftcrystal_properties_wf python ${CMAKE_CURRENT_SOURCE_DIR}/wf.py)

  endif(do_external_program_tests)

endif(do_tests)
