set(PYVASP_PY __init__.py extract.py functional.py keywords.py)
install(FILES ${PYVASP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/gulp)

add_subdirectory(tests)
