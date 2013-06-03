install( FILES __init__.py json.py extract.py makeclass.py uirepr.py 
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/tools )

add_subdirectory(input)
add_subdirectory(tests)
