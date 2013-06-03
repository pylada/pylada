set(PYVASPEXTRACT_PY __init__.py mixin.py base.py)
install(FILES ${PYVASPEXTRACT_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vasp/extract)

add_subdirectory(tests)
