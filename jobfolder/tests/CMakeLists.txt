if(do_tests)

  add_test(jobfolder python ${CMAKE_CURRENT_SOURCE_DIR}/jobfolder.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(jobparams python ${CMAKE_CURRENT_SOURCE_DIR}/jobparams.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(collect python ${CMAKE_CURRENT_SOURCE_DIR}/collect.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(forwardingdict python ${CMAKE_CURRENT_SOURCE_DIR}/forwardingdict.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(do_tests)
