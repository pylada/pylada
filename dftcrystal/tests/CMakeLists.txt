if(do_tests)
  # All will require external program to run tests.
  # This is because even the most basic operation on crystal structures runs the
  # program in the background
  if(do_external_program_tests)

    add_test(dftcrystal_rawgeom python ${CMAKE_CURRENT_SOURCE_DIR}/rawgeometry.py)
    add_test(dftcrystal_rawbasis python ${CMAKE_CURRENT_SOURCE_DIR}/rawbasis.py)
    add_test(dftcrystal_basis python ${CMAKE_CURRENT_SOURCE_DIR}/basis.py)
    add_test(dftcrystal_optgeom python ${CMAKE_CURRENT_SOURCE_DIR}/optgeom.py)
    add_test(dftcrystal_hamiltonian python ${CMAKE_CURRENT_SOURCE_DIR}/hamiltonian.py)
    add_test(dftcrystal_electronic python ${CMAKE_CURRENT_SOURCE_DIR}/electronic.py)
    add_test(dftcrystal_geometry python ${CMAKE_CURRENT_SOURCE_DIR}/geometry.py)
    add_test(dftcrystal_inputgeom python ${CMAKE_CURRENT_SOURCE_DIR}/inputgeom.py)
    add_test(dftcrystal_pickles python ${CMAKE_CURRENT_SOURCE_DIR}/pickles.py)
    add_test(dftcrystal_external python ${CMAKE_CURRENT_SOURCE_DIR}/external.py)
    add_test(dftcrystal_read_molecule python ${CMAKE_CURRENT_SOURCE_DIR}/read_molecule.py)
    add_test(dftcrystal_input python ${CMAKE_CURRENT_SOURCE_DIR}/input.py)
    add_test(dftcrystal_input_structure_regression python ${CMAKE_CURRENT_SOURCE_DIR}/input_structure_regression.py)
    add_test(dftcrystal_input_structure_regression python ${CMAKE_CURRENT_SOURCE_DIR}/input_structure_regression1.py)
    add_test(dftcrystal_gulp python ${CMAKE_CURRENT_SOURCE_DIR}/gulp.py)
    add_test(dftcrystal_crystal_regression python ${CMAKE_CURRENT_SOURCE_DIR}/crystal_regression.py)
    add_test(dftcrystal_extract_structure python ${CMAKE_CURRENT_SOURCE_DIR}/extract_structure.py)
    add_test(dftcrystal_read_input python ${CMAKE_CURRENT_SOURCE_DIR}/read_input.py)
    add_test(dftcrystal_run python ${CMAKE_CURRENT_SOURCE_DIR}/run.py)
    add_test(dftcrystal_runrelax python ${CMAKE_CURRENT_SOURCE_DIR}/runrelax.py)

  endif(do_external_program_tests)

endif(do_tests)
