install( FILES __init__.py block.py keywords.py tree.py listblock.py
         DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/tools/input )
add_subdirectory(tests)
