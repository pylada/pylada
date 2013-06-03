set(PYDFTCRYSTAL_PY __init__.py keywords.py extract.py emass.py) 

install(FILES ${PYDFTCRYSTAL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/dftcrystal/properties)

add_subdirectory(tests)
