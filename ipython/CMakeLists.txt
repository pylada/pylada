set(PYIPY_PY __init__.py explore.py savefolders.py goto.py listfolders.py
             showme.py export.py jmol.py manipfolders.py)# qstat.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/ipython)
add_subdirectory(launch)

add_subdirectory(tests)
