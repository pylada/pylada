#include "PyladaConfig.h"

#include <Python.h>
#include <structmember.h>
#include "include_numpy.h"

#include <errors/exceptions.h>
#include <math/math.h>

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif
#define PYLADA_PYTHON_MODULE 0
#include "python.h"

namespace Pylada
{
  namespace python
  {
    namespace 
    {
#     include "object.cc"
#     include "quantity.cc"
    }
  }
}



PyMODINIT_FUNC initcppwrappers(void) 
{
  using namespace Pylada::python;
  using namespace Pylada;
  static void *api_capsule[PYLADA_SLOT(python)];
  static PyMethodDef methods_table[] = { {NULL, NULL, 0, NULL} };
  PyObject *c_api_object;

  char const doc[] =  "Wrapper around basic C++ helper functions.\n\n"
                      "This module only contains a capsule for cpp functions.\n";
  PyObject* module = Py_InitModule3("cppwrappers", methods_table, doc);
  if(not module) return;

  import_array(); // needed for NumPy 

  /* Initialize the C API pointer array */
# undef PYLADA_PYTHON_PYTHON_H
# include "python.h"

  /* Create a Capsule containing the API pointer array's address */
  # ifdef PYLADA_PYTHONTWOSIX
      c_api_object = PyCObject_FromVoidPtr((void *)api_capsule, NULL);
  # else
      static const char name[] = "pylada.cppwrappers._C_API";
      c_api_object = PyCapsule_New((void *)api_capsule, name, NULL);
  # endif
  if (c_api_object != NULL) PyModule_AddObject(module, "_C_API", c_api_object);
}
