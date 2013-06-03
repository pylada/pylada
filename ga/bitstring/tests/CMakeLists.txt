if(do_tests)
  add_test(ga_bitstring python ${CMAKE_CURRENT_SOURCE_DIR}/test.py
                            ${CMAKE_CURRENT_BINARY_DIR}/../../../process/tests/pifunc)
endif(do_tests)

