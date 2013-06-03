#include "PyladaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL pylada_ewald_ARRAY_API
#include <python/include_numpy.h>

#include <algorithm>

#include <errors/exceptions.h>
#include <python/python.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#include "ewald.h"

namespace Pylada
{
  namespace pcm
  {
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"ewald", (PyCFunction)ewald, METH_KEYWORDS,
         "Performs an ewald summation." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}

PyMODINIT_FUNC initcppwrappers(void) 
{
  import_array(); // needed for NumPy 

  char const doc[] =  "Wrapper around C++/fortan point-ion models methods.";
  PyObject* module = Py_InitModule3("cppwrappers", Pylada::pcm::methods_table, doc);
  if(not module) return;
  import_array();
  if(not Pylada::python::import()) return;
}
