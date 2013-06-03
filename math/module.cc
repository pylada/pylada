#include "PyladaConfig.h"

#include <Python.h>
#include <python/include_numpy.h>

#include <python/python.h>

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif
#define PYLADA_MATH_MODULE 0
#include "math.h"

namespace Pylada
{
  namespace math
  {
    namespace details
    {
      //! \brief function which needs to be compiled without optimization.
      //! \throw error::singular_matrix if the matrix is singular.
      bool no_opt_change_test(types::t_real _new, types::t_real _last);
    }

    namespace 
    {
#     include "gruber.cc"
#     include "smith_normal_form.cc"
#     include "methods.cc"
    }
  }
}



PyMODINIT_FUNC initmath(void) 
{
  using namespace Pylada::math;
  using namespace Pylada;
  static void *api_capsule[PYLADA_SLOT(math)];
  PyObject *c_api_object;

  char const doc[] =  "Wrapper around basic C++ helper functions.\n\n"
                      "This module only contains a capsule for cpp functions.\n";
  PyObject* module = Py_InitModule3("math", methods_table, doc);
  if(not module) return;
  if(not Pylada::python::import()) return;
  import_array(); // imported by python

  /* Initialize the C API pointer array */
# undef PYLADA_MATH_MODULE_H
# include "math.h"

  /* Create a Capsule containing the API pointer array's address */
# ifdef PYLADA_PYTHONTWOSIX
    c_api_object = PyCObject_FromVoidPtr((void *)api_capsule, NULL);
# else
    static const char name[] = "pylada.math._C_API";
    c_api_object = PyCapsule_New((void *)api_capsule, name, NULL);
# endif
  if (c_api_object != NULL) PyModule_AddObject(module, "_C_API", c_api_object);
}
