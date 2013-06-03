set(PYPROCEED_PY call.py  __init__.py  iterator.py  jobfolder.py 
                 process.py  program.py mpi.py pool.py dummy.py ) 
install(FILES ${PYPROCEED_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/process)

add_subdirectory(tests)
