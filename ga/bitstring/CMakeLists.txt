set(PYGA_PY __init__.py operators.py)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ga/bitstring)

add_subdirectory(tests)
