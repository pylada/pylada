#include "PyladaConfig.h"

#include <Python.h>
#include <python/include_numpy.h>
#include "../crystal.h"

using namespace Pylada::crystal;
PyObject* get_static_object(PyObject* _module, PyObject*)
{ 
  Pylada::python::Object module = PyImport_ImportModule("_atom_self");
  if(not module) return NULL;
  PyObject *result = PyObject_GetAttrString(module.borrowed(), "_atom");
  return result;
}
PyObject* set_static_object(PyObject* _module, PyObject *_object)
{
  if(not check_atom(_object))
  {
    PYLADA_PYERROR(TypeError, "Wrong type.");
    return NULL;
  }
  Pylada::python::Object module = PyImport_ImportModule("_atom_self");
  if(not module) return NULL;
  if( PyObject_SetAttrString(module.borrowed(), "_atom", _object) < 0) return NULL;
  Py_RETURN_NONE;
}

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#ifdef PYLADA_DECLARE
#  error PYLADA_DECLARE already defined.
#endif
#define PYLADA_DECLARE(name, args) {#name, (PyCFunction)name, METH_ ## args, ""} 

static PyMethodDef methods[] = { 
  PYLADA_DECLARE(get_static_object, NOARGS),
  PYLADA_DECLARE(set_static_object, O),
  {NULL},
};

#undef PYLADA_DECLARE

PyMODINIT_FUNC init_atom_self(void) 
{
  PyObject* module = Py_InitModule("_atom_self", methods);
  if(not module) return; 
  if(not Pylada::python::import()) return;
  if(not Pylada::crystal::import()) return;
  Atom satom; 
  PyModule_AddObject(module, "_atom", (PyObject *)satom.new_ref());
  satom.release();
}
