set(PYVASP_PY files.py __init__.py specie.py functional.py relax.py emass.py
              keywords.py )

install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vasp)

# add_subdirectory(incar)
add_subdirectory(extract)
add_subdirectory(nlep)

add_subdirectory(tests)
