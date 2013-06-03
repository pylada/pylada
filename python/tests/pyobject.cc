#include "PyladaConfig.h"

#include <Python.h>
#include "../include_numpy.h"

#include "../python.h"


using namespace Pylada::python;
PyObject* represent(PyObject *_module, PyObject *_in)
{ 
  Object o = Object::acquire(_in);
  std::ostringstream sstr;
  try { sstr << o; }
  catch(std::exception &_e)
  {
    PYLADA_PYERROR(internal, "caught error");
    return NULL;
  }
  return PyString_FromString(sstr.str().c_str());
}

PyObject* add_attribute(PyObject *_module, PyObject* _args)
{
  Object o = Object::acquire(PyTuple_GET_ITEM(_args, 0));
  PyObject_SetAttr(o.borrowed(), PyTuple_GET_ITEM(_args, 1), PyTuple_GET_ITEM(_args, 2));
  Py_RETURN_NONE;
}
PyObject* callme(PyObject *_module, PyObject* _in)
{
  Object o = Object::acquire(_in);
  return PyObject_CallObject(o.borrowed(), NULL);
}

PyObject* equality(PyObject* _module, PyObject* _args)
{
  Object a = Object::acquire(PyTuple_GET_ITEM(_args, 0));
  Object b = Object::acquire(PyTuple_GET_ITEM(_args, 1));
  if(a == b) Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#ifdef PYLADA_DECLARE
#  error PYLADA_DECLARE already defined.
#endif
#define PYLADA_DECLARE(name, args) {#name, (PyCFunction)name, METH_ ## args, ""} 

static PyMethodDef methods[] = { 
  PYLADA_DECLARE(represent, O),
  PYLADA_DECLARE(add_attribute, VARARGS),
  PYLADA_DECLARE(callme, O),
  PYLADA_DECLARE(equality, VARARGS),
  {NULL},
};

#undef PYLADA_DECLARE

PyMODINIT_FUNC init_pyobject(void) 
{
  PyObject* module = Py_InitModule("_pyobject", methods);
  if(not module) return;
  if(not Pylada::python::import()) return;
}
