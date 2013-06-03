#include "PyladaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL pylada_vff_ARRAY_API
#define NO_IMPORT_ARRAY
#include <python/include_numpy.h>


#define PYLADA_NO_IMPORT
#include <errors/exceptions.h>
#include <crystal/crystal.h>

#include "../node/pybase.h"
#include "pybase.h"

namespace Pylada
{
  namespace vff
  {
    //! Creates a new edge.
    EdgeData* PyEdge_New()
    {
      EdgeData* result = (EdgeData*) edge_type()->tp_alloc(edge_type(), 0);
      if(not result) return NULL;
      new(&result->translation) math::rVector3d(0,0,0);
      result->do_translate = false;
      result->a = NULL;
      result->b = NULL;
      return result;
    }

    //! Creates a new structure with a given type.
    static EdgeData* new_with_args(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
      { return PyEdge_New(); }

    static int traverse(EdgeData *self, visitproc visit, void *arg)
    {
      if(self->a != self->b) Py_VISIT(self->b);
      Py_VISIT(self->a);
      return 0;
    }
  
    static int gcclear(EdgeData *self)
    { 
      if(self->a != self->b) Py_CLEAR(self->b);
      Py_CLEAR(self->a);
      return 0;
    }

    // Function to deallocate a string atom.
    static void dealloc(EdgeData *_self)
    {
      gcclear(_self);
  
      // calls destructor explicitely.
      PyTypeObject* ob_type = _self->ob_type;
      _self->~EdgeData();

      ob_type->tp_free((PyObject*)_self);
    }

    PyObject* edge_to_tuple(EdgeData* _data)
    {
      python::Object result = PyTuple_New(3);
      if(not result) return NULL;
      PyObject* translation = python::numpy::wrap_to_numpy(_data->translation, (PyObject*)_data);
      if(not translation) return NULL;
      Py_INCREF(_data->a);
      Py_INCREF(_data->b);
      PyTuple_SET_ITEM(result.borrowed(), 0, (PyObject*)_data->a);
      PyTuple_SET_ITEM(result.borrowed(), 1, (PyObject*)_data->b);
      PyTuple_SET_ITEM(result.borrowed(), 2, translation);
      return result.release();
    }
    PyObject* edge_to_tuple(NodeData* _self, EdgeData* _data)
    {
      python::Object result = PyTuple_New(2);
      if(not result) return NULL;
      bool const isself = _self == _data->a;
      PyObject* translation;
      if(isself)
      {
        math::rVector3d const &trans = _data->translation;
        translation = python::numpy::wrap_to_numpy(trans, (PyObject*) _data);
      }
      else
      {
        typedef math::rVector3d::Scalar t_type;
        int const nptype = python::numpy::type<t_type>::value;
        npy_intp dims[1] = {3};
        translation = PyArray_SimpleNew(1, dims, nptype);
        if(not translation) return NULL;
        *(t_type*)PyArray_GETPTR1(translation, 0) = -_data->translation[0];
        *(t_type*)PyArray_GETPTR1(translation, 1) = -_data->translation[1];
        *(t_type*)PyArray_GETPTR1(translation, 2) = -_data->translation[2];
      }
      if(not translation) return NULL;
      Py_INCREF(isself? _data->b: _data->a);
      PyTuple_SET_ITEM(result.borrowed(), 0, isself? (PyObject*)_data->b: (PyObject*)_data->a);
      PyTuple_SET_ITEM(result.borrowed(), 1, translation);
      return result.release();
    }
  
    // Returns pointer to edge type.
    PyTypeObject* edge_type()
    {
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "pylada.vff.cppwrappers.Edge",       /*tp_name*/
          sizeof(EdgeData),                  /*tp_basicsize*/
          0,                                 /*tp_itemsize*/
          (destructor)dealloc,               /*tp_dealloc*/
          0,                                 /*tp_print*/
          0,                                 /*tp_getattr*/
          0,                                 /*tp_setattr*/
          0,                                 /*tp_compare*/
          0,                                 /*tp_repr*/
          0,                                 /*tp_as_number*/
          0,                                 /*tp_as_sequence*/
          0,                                 /*tp_as_mapping*/
          0,                                 /*tp_hash */
          0,                                 /*tp_call*/
          0,                                 /*tp_str*/
          0,                                 /*tp_getattro*/
          0,                                 /*tp_setattro*/
          0,                                 /*tp_as_buffer*/
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
          "Defines a edge.\n\n",             /*tp_doc*/
          (traverseproc)traverse,            /* tp_traverse */
          (inquiry)gcclear,                  /* tp_clear */
          0,		                             /* tp_richcompare */
          0,                                 /* tp_weaklistoffset */
          0,                                 /* tp_iter */
          0,		                             /* tp_iternext */
          0,                                 /* tp_methods */
          0,                                 /* tp_members */
          0,                                 /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          0,                                 /* tp_dictoffset */
          0,                                 /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)new_with_args,            /* tp_new */
      };
      return &dummy;
    }

  } // namespace vff
} // namespace Pylada
