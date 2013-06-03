if(do_tests)

  add_test(vff_node python ${CMAKE_CURRENT_SOURCE_DIR}/node.py)
  add_test(vff_tree python ${CMAKE_CURRENT_SOURCE_DIR}/tree.py)
  add_test(vff_iterator python ${CMAKE_CURRENT_SOURCE_DIR}/iterator.py)
  add_test(vff_represent python ${CMAKE_CURRENT_SOURCE_DIR}/represent.py)
  add_test(vff_ternary python ${CMAKE_CURRENT_SOURCE_DIR}/ternary.py)
  add_test(vff_gradients python ${CMAKE_CURRENT_SOURCE_DIR}/gradients.py)
  add_test(vff_gradientfunc python ${CMAKE_CURRENT_SOURCE_DIR}/gradientfunc.py)
  add_test(vff_sigegradients python ${CMAKE_CURRENT_SOURCE_DIR}/sigegradients.py)
  add_test(vff_functional python ${CMAKE_CURRENT_SOURCE_DIR}/functional.py)
  add_test(vff_epirelax python ${CMAKE_CURRENT_SOURCE_DIR}/epirelax.py)

endif(do_tests)
