if(do_tests)

  add_test(vasp_extraxt_common python ${CMAKE_CURRENT_SOURCE_DIR}/common.py ${CMAKE_CURRENT_SOURCE_DIR}/data)
  add_test(vasp_extraxt_dft python ${CMAKE_CURRENT_SOURCE_DIR}/dft.py ${CMAKE_CURRENT_SOURCE_DIR}/data)
  
endif(do_tests)

