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
       //! Adds atom to structure.
       static PyObject* link(NodeData* _self, PyObject* _args, PyObject* _kwargs);
       //! Removes and returns key from set.
       static PyObject* pop(NodeData *_self, PyObject* _index);
       //! Clears occupation set.
       static PyObject* clear(NodeData *_self);
       //! Returns size of set.
       static Py_ssize_t size(NodeData* _self) { return Py_ssize_t(_self->bonds.size()); }
       //! Retrieves item. No slice.
       static PyObject* getitem(NodeData *_self, Py_ssize_t _index);
       //! Retrieves slice.
       static PyObject* subscript(NodeData *_self, PyObject *_index);
    }
#   ifdef PYLADA_STARTTRY
#     error PYLADA_STARTTRY already defined.
#   endif
#   ifdef PYLADA_ENDTRY
#     error PYLADA_ENDTRY already defined.
#   endif
#   define PYLADA_STARTTRY try {
#   define PYLADA_ENDTRY(error) } catch(std::exception &_e) { \
       PYLADA_PYERROR_FORMAT( internal, \
                            "C++ thrown exception caught in " #error ":\n%.200s", \
                            _e.what() ); \
       return NULL; } 

    // Adds atom to structure.
    PyObject* link(NodeData* _self, PyObject* _args, PyObject* _kwargs)
    {
      static char *kwlist[] = { const_cast<char*>("endpoint"), 
                                const_cast<char*>("translation"),
                                NULL };
      PyObject *nodeB = NULL;
      PyObject *_pos = NULL;
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "O|O:Node.append", kwlist,
                                          &nodeB, &_pos ) )
        return NULL;
      if(not PyNodeData_Check(nodeB))
      {
        PYLADA_PYERROR(TypeError, "First argument to link should be a Node.");
        return NULL;
      }
      PYLADA_STARTTRY
        size_t const n = _self->bonds.size();
        math::rVector3d pos(0,0,0);
        if(_pos) 
          { if(not python::numpy::convert_to_vector(_pos, pos)) return NULL; }
        if(not PyNode_AddEdge(_self, (NodeData*)nodeB, pos)) return NULL;
        if(n != _self->bonds.size()) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
      PYLADA_ENDTRY(Node.link)
    }
    // Removes and returns atom from set.
    PyObject* pop(NodeData *_self, PyObject* _index)
    {
      PYLADA_STARTTRY
        long index = PyInt_AsLong(_index);
        if(index == -1 and PyErr_Occurred()) return NULL;
        size_t const n = _self->bonds.size();
        if(index < 0) index += n;
        if(index < 0 or index >= n) 
        {
          PYLADA_PYERROR(IndexError, "Index out-of-range in Node.pop.");
          return NULL;
        }

        // edge variable takes ownership.
        python::Object edge((PyObject*)_self->bonds[index]);
        EdgeData* const py_edge = (EdgeData*) edge.borrowed();
        _self->bonds.erase(_self->bonds.begin()+index);

        //! Remove bond from other atom.
        if(py_edge->a != py_edge->b)
        {
          NodeData * other = py_edge->a != _self ? py_edge->a: py_edge->b;
          std::vector<EdgeData*> :: iterator i_first = other->bonds.begin();
          std::vector<EdgeData*> :: iterator const i_end = other->bonds.end();
          for(; i_first != i_end; ++i_first)
            if((*i_first) == py_edge) 
            {
              other->bonds.erase(i_first);
              Py_DECREF(py_edge);
              break;
            }
          if(i_first == i_end)
          {
            PYLADA_PYERROR(internal, "Could not find equivalent edge in other endpoint.");
            return NULL;
          }
        }

        return edge_to_tuple(_self, py_edge);
      PYLADA_ENDTRY(Node.pop)
    }
    // Clears occupation set.
    PyObject* clear(NodeData *_self)
    {
      // hold a reference to make sure self is not deleted until the very end.
      python::Object const hold = python::Object::acquire((PyObject*)_self); 
      std::vector<EdgeData*> :: const_iterator i_bond = _self->bonds.begin();
      std::vector<EdgeData*> :: const_iterator i_bond_end = _self->bonds.end();
      for(; i_bond != i_bond_end; ++i_bond)
      {
        if((*i_bond)->a == (*i_bond)->b) continue;
        NodeData * const other = (*i_bond)->a != _self ? (*i_bond)->a: (*i_bond)->b;
        std::vector<EdgeData*> :: iterator i_first = other->bonds.begin();
        std::vector<EdgeData*> :: iterator const i_end = other->bonds.end();
        for(; i_first != i_end; ++i_first)
          if(*i_first == *i_bond) 
          {
            other->bonds.erase(i_first);
            Py_DECREF((*i_bond));
            break;
          }
        if(i_first == i_end)
        {
          PYLADA_PYERROR(internal, "Could not find equivalent edge in other endpoint.");
          return NULL;
        }
      }
      // copy vector so we can delete bonds.
      std::vector<EdgeData*> const copy(_self->bonds);
      // erase all bonds.
      _self->bonds.clear();
      i_bond = copy.begin();
      i_bond_end = copy.end();
      for(; i_bond != i_bond_end; ++i_bond) { Py_DECREF((*i_bond)); }
      Py_RETURN_NONE;
    }
    // Retrieves item. No slice.
    PyObject* getitem(NodeData *_self, Py_ssize_t _index)
    {
      PYLADA_STARTTRY
        size_t const n = _self->bonds.size();
        if(_index < 0) _index += n;
        if(_index < 0 or _index >= n)
        {
          PYLADA_PYERROR(IndexError, "Index out of range when getting node from structure.");
          return NULL;
        }
        EdgeData* const py_edge = _self->bonds[_index];
        return edge_to_tuple(_self, py_edge);
      PYLADA_ENDTRY(Node.__getitem__)
    }
    // Retrieves slice and items.
    PyObject* subscript(NodeData *_self, PyObject *_index)
    {
      if(PyIndex_Check(_index))
      {
        Py_ssize_t i = PyNumber_AsSsize_t(_index, PyExc_IndexError);
        if (i == -1 and PyErr_Occurred()) return NULL;
        return getitem(_self, i);
      }
      else if(not PySlice_Check(_index))
      {
        PYLADA_PYERROR_FORMAT( TypeError,
                             "Structure indices must be integers or slices, not %.200s.\n",
                             _index->ob_type->tp_name );
        return NULL;
      }
      PYLADA_STARTTRY
        Py_ssize_t start, stop, step, slicelength;
        if (PySlice_GetIndicesEx((PySliceObject*)_index, _self->bonds.size(),
                                  &start, &stop, &step, &slicelength) < 0)
          return NULL;
        if(slicelength == 0) return PyTuple_New(0);
        python::Object tuple = PyTuple_New(slicelength);
        if(not tuple) return NULL;
        if(step < 0)
        {
          std::vector<EdgeData*>::const_reverse_iterator i_first = _self->bonds.rbegin()
                                                 + (_self->bonds.size() - start-1);
          for(Py_ssize_t i(0); i < slicelength; i_first -= step, ++i)
          {
            PyObject *bond = edge_to_tuple(_self, *i_first);
            if(not bond) return NULL;
            PyTuple_SET_ITEM(tuple.borrowed(), i, bond);
          }
        }
        else 
        {
          std::vector<EdgeData*>::const_iterator i_first = _self->bonds.begin() + start;
          for(Py_ssize_t i(0); i < slicelength; i_first += step, ++i)
          {
            PyObject *bond = edge_to_tuple(_self, *i_first);
            if(not bond) return NULL;
            PyTuple_SET_ITEM(tuple.borrowed(), i, bond);
          }
        }
        return tuple.release();
      PYLADA_ENDTRY(Node.__getitem__);
    }
#   undef PYLADA_STARTTRY
#   undef PYLADA_ENDTRY
  }
}
