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

namespace Pylada
{
  namespace vff
  {
    extern "C" 
    {
      //! Describes a node in a first neighbor net.
      struct BondIterator
      {
        PyObject_HEAD 
        //! Center of the bonds which are being iterated.
        NodeData* parent;
        //! First and last iterator.
        std::vector<EdgeData*> :: iterator i_first, i_end;
        //! Whether the first next was called. 
        bool is_first; 
        //! Size on creation.
        std::vector<EdgeData*> :: size_type N;
      };
      //! Returns Self.
      static PyObject* bond_iterator_iter(PyObject* _self) { Py_INCREF(_self); return _self; }
      //! Returns next object.
      static PyObject* bond_iterator_next(BondIterator* _self);
      //! Function to deallocate a string atom.
      static void bond_iterator_dealloc(BondIterator *_self);
      //! Creates iterator.
      static PyObject* bond_iterator_create(NodeData* _self);
      //! Returns next object without double counting.
      static PyObject* dcbond_iterator_next(BondIterator* _self);
      
      //! Describes a node in a first neighbor net.
      struct AngleIterator
      {
        PyObject_HEAD 
        //! Center of the bonds which are being iterated.
        NodeData* parent;
        //! Tuples to be yielded.
        std::vector<PyObject*> container;
        //! First and last iterator.
        std::vector<PyObject*> :: iterator i_first, i_second, i_end;
      };
      //! Returns next object.
      static PyObject* angle_iterator_next(AngleIterator* _self);
      //! Function to deallocate a string atom.
      static void angle_iterator_dealloc(AngleIterator *_self);
      //! Creates iterator.
      static PyObject* angle_iterator_create(NodeData* _self);
    }

    // Returns next object.
    PyObject* bond_iterator_next(BondIterator* _self)
    {
      if(_self->N != _self->parent->bonds.size())
      {
        PYLADA_PYERROR(IndexError, "Bonds changed mid-iteration.");
        return NULL;
      }
      if(_self->i_first != _self->i_end) 
      {
        if(_self->is_first) _self->is_first = false; else ++_self->i_first;
        if(_self->i_first != _self->i_end)
          return edge_to_tuple(_self->parent, *_self->i_first);
      }
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    //! Function to deallocate a string atom.
    static void bond_iterator_dealloc(BondIterator *_self)
    {
      PyObject* dummy = (PyObject*)_self->parent;
      _self->parent = NULL;
      Py_XDECREF(dummy);
    }
    //! Creates iterator.
    template<PyTypeObject*(*itertype)()>
      static PyObject* bond_iterator_create(NodeData* _in)
      {
        BondIterator *result = PyObject_New(BondIterator, itertype());
        if(result == NULL) return NULL;
        result->parent = _in;
        Py_INCREF(result->parent);
        new(&result->i_first) std::vector<EdgeData*>::iterator(_in->bonds.begin());
        new(&result->i_end) std::vector<EdgeData*>::iterator(_in->bonds.end());
        result->is_first = true;
        result->N = _in->bonds.size();
        return (PyObject*) result;
      }
    //! Returns self.
    static PyObject* get_self(PyObject* _self) { Py_INCREF(_self); return _self; }

    // Returns next object.
    PyObject* dcbond_iterator_next(BondIterator* _self)
    {
      if(_self->N != _self->parent->bonds.size())
      {
        PYLADA_PYERROR(IndexError, "Bonds changed mid-iteration.");
        return NULL;
      }
      if(_self->i_first != _self->i_end) 
      {
        if(_self->is_first) _self->is_first = false; else ++_self->i_first;
        while(_self->i_first != _self->i_end)
        {
          if((*_self->i_first)->a == _self->parent)
            return edge_to_tuple(_self->parent, *_self->i_first);
          ++_self->i_first;
        }
      }
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }

    // Returns next object.
    PyObject* angle_iterator_next(AngleIterator* _self)
    {
      if(_self->i_first == _self->i_end)
      {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
      }
      if(++_self->i_second == _self->i_end) 
      {
        _self->i_second = ++_self->i_first;
        return angle_iterator_next(_self);
      }
      
      PyObject *result = PyTuple_New(2);
      Py_INCREF(*_self->i_first);
      PyTuple_SET_ITEM(result, 0, *_self->i_first);
      Py_INCREF(*_self->i_second);
      PyTuple_SET_ITEM(result, 1, *_self->i_second);
      return result;
    }
    //! Function to deallocate a string atom.
    static void angle_iterator_dealloc(AngleIterator *_self)
    {
      PyObject* dummy = (PyObject*)_self->parent;
      _self->parent = NULL;
      Py_XDECREF(dummy);
      std::vector<PyObject*>::iterator i_first = _self->container.begin();
      std::vector<PyObject*>::iterator const i_end = _self->container.end();
      for(; i_first != i_end; ++i_first)
      {
        dummy = *i_first;
        *i_first = NULL;
        Py_XDECREF(dummy);
      }
      _self->container.clear();
      Py_XDECREF(dummy);
    }
    //! Creates iterator.
    static PyObject* angle_iterator_create(NodeData* _in)
    {
      AngleIterator *result = PyObject_New(AngleIterator, angleiterator_type());
      if(result == NULL) return NULL;
      result->parent = _in;
      Py_INCREF(result->parent);
      new(&result->container) std::vector<PyObject*>;
      std::vector<EdgeData*>::iterator i_first = _in->bonds.begin();
      std::vector<EdgeData*>::iterator const i_end = _in->bonds.end();
      for(; i_first != i_end; ++i_first)
      {
        PyObject *dummy = edge_to_tuple(_in, *i_first);
        if(dummy == NULL) return NULL;
        result->container.push_back(dummy);
      }
      new(&result->i_first) std::vector<PyObject*>::iterator(result->container.begin());
      new(&result->i_second) std::vector<PyObject*>::iterator(result->i_first);
      new(&result->i_end) std::vector<PyObject*>::iterator(result->container.end());
      return (PyObject*) result;
    }
  }
}
