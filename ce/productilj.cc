/******************************
   This file is part of PyLaDa.

   Copyright (C) 2013 National Renewable Energy Lab
  
   PyLaDa is a high throughput computational platform for Physics. It aims to make it easier to submit
   large numbers of jobs on supercomputers. It provides a python interface to physical input, such as
   crystal structures, as well as to a number of DFT (VASP, CRYSTAL) and atomic potential programs. It
   is able to organise and launch computational jobs on PBS and SLURM.
  
   PyLaDa is free software: you can redistribute it and/or modify it under the terms of the GNU General
   Public License as published by the Free Software Foundation, either version 3 of the License, or (at
   your option) any later version.
  
   PyLaDa is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
   Public License for more details.
  
   You should have received a copy of the GNU General Public License along with PyLaDa.  If not, see
   <http://www.gnu.org/licenses/>.
******************************/

#include "PyladaConfig.h"

#include <Python.h>
#include <errors/exceptions.h>

#include "productilj.h"

namespace Pylada
{
  namespace ce
  {
    //! Creates a new ndimiterator.
    ProductILJIterator* PyProductILJIterator_New()
    {
      ProductILJIterator* result = (ProductILJIterator*) productiljiterator_type()
        ->tp_alloc(productiljiterator_type(), 0);
      if(not result) return NULL;
      result->sequence = NULL;
      new(&result->counter) std::vector<Py_ssize_t>;
      return result;
    }
    //! Creates a new ndimiterator with a given type.
    ProductILJIterator* PyProductILJIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      ProductILJIterator* result = (ProductILJIterator*)_type->tp_alloc(_type, 0);
      if(not result) return NULL;
      result->sequence = NULL;
      new(&result->counter) std::vector<Py_ssize_t>;
      return result;
    }

    // Creates a new ndimiterator with a given type, also calling initialization.
    ProductILJIterator* PyProductILJIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    {
      ProductILJIterator* result = PyProductILJIterator_NewWithArgs(_type, _args, _kwargs);
      if(result == NULL) return NULL;
      if(_type->tp_init((PyObject*)result, _args, _kwargs) < 0) {Py_DECREF(result); return NULL; }
      return result;
    }

    // Function to deallocate a string atom.
    static void dealloc(ProductILJIterator *_self)
    {
      if(_self->sequence != NULL)
      { 
        PyObject *dummy = _self->sequence;
        _self->sequence = NULL;
        Py_DECREF(dummy);
      }

      // calls destructor explicitely.
      PyTypeObject* ob_type = _self->ob_type;
      _self->~ProductILJIterator();

      ob_type->tp_free((PyObject*)_self);
    }
  
    // Function to initialize an atom.
    static int init(ProductILJIterator* _self, PyObject* _args, PyObject *_kwargs)
    {
      unsigned long i = 0;
      PyObject *object = NULL;
      if(not PyArg_ParseTuple(_args, "Ok:ProductILJ", &object, &i))
         return -1;
      if(_self->sequence != NULL)
      {
        PyObject *dummy = _self->sequence;
        _self->sequence = NULL;
        Py_DECREF(dummy);
      }
      if(PySequence_Check(object))
      {
        _self->sequence = PySequence_List(object);
        if(not _self->sequence) return -1;
      }
      else if(PyIter_Check(object))
      {
        _self->sequence = PyList_New(0);
        if(not _self->sequence) return -1;
        while(PyObject *item = PyIter_Next(object))
          PyList_Append(_self->sequence, item);
        if(PyErr_Occurred() != NULL) return -1;
      }
      else 
      {
        PYLADA_PYERROR(TypeError, "First argument to ProductILJ is not iterable.");
        return -1;
      }
      _self->N = PyList_Size(_self->sequence);
      _self->counter.clear();
      _self->counter.resize(i);
      std::vector<Py_ssize_t>::iterator i_first = _self->counter.begin();
      std::vector<Py_ssize_t>::iterator const i_end = _self->counter.end();
      for(Py_ssize_t j(0); i_first != i_end; ++i_first, ++j) *i_first = j;
      _self->is_first = true;
      return 0;
    }
  
    static PyObject* self(PyObject* _self)
    {
      Py_INCREF(_self);
      return _self;
    }

    static Py_ssize_t increment(
          std::vector<Py_ssize_t>::reverse_iterator const &_current, 
          std::vector<Py_ssize_t>::reverse_iterator const &_end,
          Py_ssize_t const _max )
    {
      if(++(*_current) >= _max)
      {
        std::vector<Py_ssize_t>::reverse_iterator i_next = _current;
        if(++i_next != _end)
          *_current = 1+increment(i_next, _end, _max-1);
      }
      return *_current;
    }


    static PyObject* next(ProductILJIterator* _self)
    {
      if(_self->counter.size() == 0 or _self->sequence == NULL) 
      {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
      }
      if(_self->is_first) _self->is_first = false;
      else increment(_self->counter.rbegin(), _self->counter.rend(), _self->N);
      if( _self->counter.back() >= _self->N)
      {
        PyObject* dummy = _self->sequence;
        _self->sequence = NULL;
        Py_DECREF(dummy);
        _self->counter.clear();
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
      }
      PyObject *tuple = PyTuple_New(_self->counter.size());
      if(not tuple) return NULL;
      std::vector<Py_ssize_t>::const_iterator i_first = _self->counter.begin();
      std::vector<Py_ssize_t>::const_iterator const i_end = _self->counter.end();
      for(Py_ssize_t i(0); i_first != i_end; ++i_first, ++i)
      {
        PyObject * const item = PyList_GET_ITEM(_self->sequence, *i_first);
        if(not item)
        {
          Py_DECREF(tuple);
          return NULL;
        }
        Py_INCREF(item);
        PyTuple_SET_ITEM(tuple, i, item);
      }
      return tuple;
    }


    // Returns pointer to ndimiterator type.
    PyTypeObject* productiljiterator_type()
    {
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "pylada.ce.cppwrappers.ProductILJ",  /*tp_name*/
          sizeof(ProductILJIterator),        /*tp_basicsize*/
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
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /*tp_flags*/
          "Defines a product iterator.\n\n"
          "Takes a sequence/iterator as input and an integer ``n``. It yields\n"
          "tuples with items :py:math:`i_0, i_1, ..i_n` from the sequence, such\n"
          "that :py:math:`i_0<i_1<...<i_n`.\n",
          0,                                 /* tp_traverse */
          0,                                 /* tp_clear */
          0,		                             /* tp_richcompare */
          0,                                 /* tp_weaklistoffset */
          (getiterfunc)self,                 /* tp_iter */
          (iternextfunc)next,                /* tp_iternext */
          0,                                 /* tp_methods */
          0,                                 /* tp_members */
          0,                                 /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          0,                                 /* tp_dictoffset */
          (initproc)init,                    /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)PyProductILJIterator_NewWithArgs,  /* tp_new */
      };
      return &dummy;
    }

  } // namespace Crystal
} // namespace Pylada
