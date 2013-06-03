set(PYGA_PY operators.py objective.py initialization.py __init__.py functional.py)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ga/xgsgo)

add_subdirectory(tests)
