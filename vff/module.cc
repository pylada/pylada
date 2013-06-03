#include "PyladaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL pylada_vff_ARRAY_API
#include <python/include_numpy.h>

#include <algorithm>

#include <errors/exceptions.h>
#include <crystal/crystal.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#include "node/pybase.h"
#include "edge/pybase.h"
namespace Pylada
{
  namespace vff
  {
#   include "vff/zb.cc"
  }
}

namespace Pylada
{
  namespace vff
  {
    //! Methods table for vff module.
    static PyMethodDef methods_table[] = {
        { "_zb_energy", (PyCFunction) zincblende::energy, METH_VARARGS,
          "Computes energy for zinc-blende\n\n"
          "This is a very specialized function to compute energies for\n"
          "zinc-blende. It expects 6 bond-parameters and 7 angle-parameters in\n"
          "numpy arrays of type 'float64'. Anything else may have unpredictable\n"
          "results." },
        { "_zb_jacobian", (PyCFunction) zincblende::jacobian, METH_VARARGS,
          "Computes jacobian for zinc-blende\n\n"
          "This is a very specialized function to compute jacobian for\n"
          "zinc-blende. It expects 6 bond-parameters and 7 angle-parameters in\n"
          "numpy arrays of type 'float64'. Anything else may have unpredictable\n"
          "results." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}
PyMODINIT_FUNC initcppwrappers(void) 
{
  char const doc[] =  "Wrapper around C++ vff class and affiliates.";
  PyObject* module = Py_InitModule3("cppwrappers", Pylada::vff::methods_table, doc);
  if(not module) return;
  import_array(); // needed for NumPy 
  if(not Pylada::python::import()) return;
  if(not Pylada::math::import()) return;
  if(not Pylada::crystal::import()) return;

  if (PyType_Ready(Pylada::vff::node_type()) < 0) return;
  if (PyType_Ready(Pylada::vff::edge_type()) < 0) return;
  if (PyType_Ready(Pylada::vff::bonditerator_type()) < 0) return;
  if (PyType_Ready(Pylada::vff::dcbonditerator_type()) < 0) return;
  if (PyType_Ready(Pylada::vff::angleiterator_type()) < 0) return;

  Py_INCREF(Pylada::vff::node_type());
  Py_INCREF(Pylada::vff::edge_type());
  Py_INCREF(Pylada::vff::bonditerator_type());
  Py_INCREF(Pylada::vff::dcbonditerator_type());
  Py_INCREF(Pylada::vff::angleiterator_type());


  PyModule_AddObject(module, "Node", (PyObject *)Pylada::vff::node_type());
  PyModule_AddObject(module, "Edge", (PyObject *)Pylada::vff::edge_type());
  PyModule_AddObject(module, "BondIterator", (PyObject *)Pylada::vff::bonditerator_type());
  PyModule_AddObject(module, "ScBondIterator", (PyObject *)Pylada::vff::dcbonditerator_type());
  PyModule_AddObject(module, "AngleIterator", (PyObject *)Pylada::vff::angleiterator_type());
}
