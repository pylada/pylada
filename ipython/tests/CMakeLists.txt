if(do_tests)

  add_test(ipython_explore ipython ${CMAKE_CURRENT_SOURCE_DIR}/explore.py ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(ipython_goto ipython ${CMAKE_CURRENT_SOURCE_DIR}/goto.py ${CMAKE_CURRENT_SOURCE_DIR})

endif(do_tests)
