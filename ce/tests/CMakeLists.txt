if(do_tests)

  add_test(ce_cluster python ${CMAKE_CURRENT_SOURCE_DIR}/cluster.py ${CMAKE_CURRENT_BINARY_DIR})
  add_test(ce_cpp python ${CMAKE_CURRENT_SOURCE_DIR}/cpp.py ${CMAKE_CURRENT_BINARY_DIR})
# add_test(ce_factory python ${CMAKE_CURRENT_SOURCE_DIR}/factory.py ${CMAKE_CURRENT_BINARY_DIR})
  
endif(do_tests)
