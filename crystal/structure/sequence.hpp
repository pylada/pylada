//! Adds atom to structure.
PyObject* structure_append(PyStructureObject* _self, PyObject* _atom);
//! Extends atoms with other list of atoms.
PyObject* structure_extend(PyStructureObject* _self, PyObject *_b);
//! Removes and returns key from set.
PyObject* structure_pop(PyStructureObject *_self, PyObject* _index);
//! Insert atom at given position
PyObject* structure_insert(PyStructureObject* _self, PyObject* _tuple);
//! Clears occupation set.
PyObject* structure_clear(PyStructureObject *_self) { _self->atoms.clear(); Py_RETURN_NONE;}
//! Returns size of set.
Py_ssize_t structure_size(PyStructureObject* _self) { return Py_ssize_t(_self->atoms.size()); }
//! Retrieves item. No slice.
PyObject* structure_getitem(PyStructureObject *_self, Py_ssize_t _index);
//! Retrieves slice and items.
PyObject* structure_subscript(PyStructureObject *_self, PyObject* _index);
//! Sets/deletes item. No slice.
int structure_setitem(PyStructureObject *_self, Py_ssize_t _index, PyObject *_replace);
//! Sets slices and items.
int structure_ass_subscript(PyStructureObject *_self, PyObject* _index, PyObject *_value);
//! Sets slices and items. Normal(not special) method.
PyObject* structure_setitemnormal(PyStructureObject *_self, PyObject* _tuple);

#ifdef LADA_STARTTRY
#  error LADA_STARTTRY already defined.
#endif
#ifdef LADA_ENDTRY
#  error LADA_ENDTRY already defined.
#endif
#define LADA_STARTTRY try {
#define LADA_ENDTRY(error) } catch(std::exception &_e) { \
    LADA_PYERROR_FORMAT( internal, \
                         "C++ thrown exception caught in " #error ":\n%.200s", \
                         _e.what() ); \
    return NULL; } 
      
// Adds atom to structure.
PyObject* structure_append(PyStructureObject* _self, PyObject* _atom)
{
  LADA_STARTTRY
    if(not Atom::check(_atom))
    {
      LADA_PYERROR_FORMAT( TypeError, 
                           "Only Atom can be appended to structure, not %.200s.",
                           _atom->ob_type->tp_name );
      return NULL;
    }
    _self->atoms.push_back(Atom::acquire_(_atom)); 
    Py_RETURN_NONE;
  LADA_ENDTRY(Structure.append)
}
//! Extends atoms with other list of atoms.
PyObject* structure_extend(PyStructureObject* _self, PyObject* _b)
{
  LADA_STARTTRY
    if(check_structure(_b))
      std::copy( ((PyStructureObject*)_b)->atoms.begin(), ((PyStructureObject*)_b)->atoms.end(),
                 std::back_inserter(_self->atoms) );
    else if(python::Object iterator = PyObject_GetIter(_b))
      for( python::Object item = PyIter_Next(iterator.borrowed());
           item.is_valid(); 
           item = PyIter_Next(iterator.borrowed()) )
      {
        if(not Atom::check(item.borrowed()))
        {
          LADA_PYERROR_FORMAT( TypeError, 
                               "Only Atom and subtypes can be assigned to structures, not %200s.",
                               item.borrowed()->ob_type->tp_name );
          return NULL;
        }
        _self->atoms.push_back(Atom::acquire_(item.borrowed()));
      }
    else 
    {
      LADA_PYERROR_FORMAT( TypeError, 
                           "Cannot append instance of %.200s to Structure",
                           _b->ob_type->tp_name );
      return NULL;
    }
    Py_RETURN_NONE;
  LADA_ENDTRY(Structure.extend)
}  
// Removes and returns atom from set.
PyObject* structure_pop(PyStructureObject *_self, PyObject* _index)
{
  LADA_STARTTRY
    long index = PyInt_AsLong(_index);
    if(index == -1 and PyErr_Occurred()) return NULL;
    if(index < 0) index += _self->atoms.size();
    if(index < 0 or index >= long(_self->atoms.size()) ) 
    {
      LADA_PYERROR(IndexError, "Index out-of-range in Structure.pop.");
      return NULL;
    }
    Atom result = _self->atoms[index];
    _self->atoms.erase(_self->atoms.begin()+index);
    return result.release();
  LADA_ENDTRY(Structure.pop)
}
// Insert atom at given position
PyObject* structure_insert(PyStructureObject* _self, PyObject* _tuple)
{
  Py_ssize_t index;
  PyObject *atom;
  if (!PyArg_ParseTuple(_tuple, "nO:insert", &index, &atom)) return NULL;
  LADA_STARTTRY
    if(not Atom::check(atom))
    {
      LADA_PYERROR_FORMAT( TypeError, 
                           "Can only Atom and subtype in Structure, not %.200s",
                           atom->ob_type->tp_name );
      return NULL;
    }
    if(index < 0) index += _self->atoms.size();
    if(index < 0 or index > long(_self->atoms.size()) ) 
    {
      LADA_PYERROR(IndexError, "Index out-of-range in Structure.pop.");
      return NULL;
    }
    if(index == long(_self->atoms.size())) _self->atoms.push_back(Atom::acquire_(atom));
    else _self->atoms.insert(_self->atoms.begin()+index, Atom::acquire_(atom));
    Py_RETURN_NONE;
  LADA_ENDTRY(Structure.insert)
}
// Retrieves item. No slice.
PyObject* structure_getitem(PyStructureObject *_self, Py_ssize_t _index)
{
  LADA_STARTTRY
    if(_index < 0) _index += _self->atoms.size();
    if(_index < 0 or _index >= long(_self->atoms.size()))
    {
      LADA_PYERROR(IndexError, "Index out of range when getting atom from structure.");
      return NULL;
    }
    return _self->atoms[_index].new_ref(); 
  LADA_ENDTRY(Structure.__getitem__)
}
// Retrieves slice and items.
PyObject* structure_subscript(PyStructureObject *_self, PyObject *_index)
{
  if(PyIndex_Check(_index))
  {
    Py_ssize_t i = PyNumber_AsSsize_t(_index, PyExc_IndexError);
    if (i == -1 and PyErr_Occurred()) return NULL;
    return structure_getitem(_self, i);
  }
  else if(not PySlice_Check(_index))
  {
    LADA_PYERROR_FORMAT( TypeError,
                         "Structure indices must be integers or slices, not %.200s.\n",
                         _index->ob_type->tp_name );
    return NULL;
  }
  LADA_STARTTRY
    Py_ssize_t start, stop, step, slicelength;
    if (PySlice_GetIndicesEx((PySliceObject*)_index, _self->atoms.size(),
                              &start, &stop, &step, &slicelength) < 0)
      return NULL;
    if(slicelength == 0) return PyTuple_New(0);
    python::Object tuple = PyTuple_New(slicelength);
    if(not tuple) return NULL;
    if(step < 0)
    {
      std::vector<Atom>::const_reverse_iterator i_first = _self->atoms.rbegin()
                                             + (_self->atoms.size() - start-1);
      for(Py_ssize_t i(0); i < slicelength; i_first -= step, ++i)
        PyTuple_SET_ITEM(tuple.borrowed(), i, i_first->new_ref());
    }
    else 
    {
      std::vector<Atom>::const_iterator i_first = _self->atoms.begin() + start;
      for(Py_ssize_t i(0); i < slicelength; i_first += step, ++i)
        PyTuple_SET_ITEM(tuple.borrowed(), i, i_first->new_ref());
    }
    return tuple.release();
  LADA_ENDTRY(Structure.__getitem__);
}
//! Sets/deletes item. No slice.
int structure_setitem(PyStructureObject *_self, Py_ssize_t _index, PyObject *_replace)
{
  try
  {
    if(_index < 0) _index += _self->atoms.size();
    if(_index < 0 or _index >= long(_self->atoms.size()))
    {
      LADA_PYERROR(IndexError, "Index out of range when setting atom in structure.");
      return -1;
    }
    // deleting.
    if(_replace == NULL)
    {
      if(_index == long(_self->atoms.size() - 1)) _self->atoms.pop_back();
      else _self->atoms.erase(_self->atoms.begin()+_index);
    }
    else if(Atom::check(_replace)) _self->atoms[_index].Object::reset(_replace);
    else 
    {
      LADA_PYERROR_FORMAT( TypeError, 
                           "Can only Atom and subtype in Structure, not %.200s",
                           _replace->ob_type->tp_name );
      return -1;
    }
    return 0;
  } 
  catch(std::exception &_e) 
  { 
     LADA_PYERROR_FORMAT( internal, 
                          "C++ thrown exception caught in Structure.__setitem__ : %200s.",
                          _e.what() );
     return -1;
  }
  return 0;
}
// Sets slices and items.
int structure_ass_subscript(PyStructureObject *_self, PyObject* _index, PyObject *_value)
{
  if(PyIndex_Check(_index))
  {
    Py_ssize_t i = PyNumber_AsSsize_t(_index, PyExc_IndexError);
    if (i == -1 and PyErr_Occurred()) return -1;
    return structure_setitem(_self, i, _value);
  }
  try
  {
    Py_ssize_t start, stop, step, slicelength;
    if (PySlice_GetIndicesEx((PySliceObject*)_index, _self->atoms.size(),
                              &start, &stop, &step, &slicelength) < 0)
        return -1;
    if(slicelength == 0) return 0; 

    if(_value == NULL) // delete.
    {
      if(step == 1) _self->atoms.erase(_self->atoms.begin()+start, _self->atoms.begin()+stop);
      else if(step > 0)
      {
        std::vector<Atom>::iterator i_first = _self->atoms.begin() + start + step * (slicelength-1);
        for(Py_ssize_t i(0); i < slicelength; ++i, i_first -= step) _self->atoms.erase(i_first);
      }
      else if(step < 0)
      {
        std::vector<Atom>::iterator i_first = _self->atoms.begin() + start;
        for(Py_ssize_t i(0); i < slicelength; ++i, i_first += step) _self->atoms.erase(i_first);
      }
    }
    else if(_self == (PyStructureObject*)_value and step == -1) // structure[::-1] = structure
      std::reverse(_self->atoms.begin()+start + (slicelength-1) * step, _self->atoms.begin()+start+1);
    else if(check_structure(_value)) // another structure.
    {
      if(long(((PyStructureObject*)_value)->atoms.size()) != slicelength)
      {
        LADA_PYERROR(ValueError, "size of right-hand side of assignement is incorrect.");
        return -1;
      }
      std::vector<Atom>::iterator i_first = _self->atoms.begin() + start;
      std::vector<Atom>::const_iterator i_other = ((PyStructureObject*)_value)->atoms.begin();
      for(Py_ssize_t i(start); i < stop; i += step, i_first += step, ++i_other)
        *i_first = *i_other;
    }
    else if(python::Object iterator = PyObject_GetIter(_value))
    {
      std::vector<Atom>::iterator i_first = _self->atoms.begin() + start;
      Py_ssize_t i(0);
      for( python::Object item = PyIter_Next(iterator.borrowed()); 
           item.is_valid();
           item.reset(PyIter_Next(iterator.borrowed())) )
      {
        if(not Atom::check(item.borrowed()))
        {
          LADA_PYERROR_FORMAT( TypeError,
                               "Only Atom and subtype can be assigned to structure, not %.200s.",
                               item.borrowed()->ob_type->tp_name );
        }
        i_first->Object::reset(item.borrowed());
        ++i;
        if(i < slicelength) i_first += step; // don't do this unless its ok, for safety.
        else break;
      }
      if(i < slicelength)
      {
        LADA_PYERROR(IndexError, "Right-hand side too small in Structure assignment.");
        return -1;
      }
      else if(PyObject* item = PyIter_Next(iterator.borrowed()))
      {
        Py_DECREF(item);
        LADA_PYERROR(IndexError, "Left-hand side too large in Structure assignment.");
        return -1;
      }
    }
    else
    {
      LADA_PYERROR_FORMAT( TypeError, 
                           "Can only assign Atoms or subtypes to Structure, not %.200s.",
                           _value->ob_type->tp_name );
      return -1;
    }
  }
  catch(std::exception &_e) 
  { 
     LADA_PYERROR_FORMAT( internal,
                          "C++ thrown exception caught in Structure.__setitem__: %.200s",
                          _e.what() );
     return -1;
  }
  return 0;
}
// Sets slices and items. Normal(not special) method.
PyObject* structure_setitemnormal(PyStructureObject *_self, PyObject* _tuple)
{
  PyObject *index, *value;
  if (not PyArg_ParseTuple(_tuple, "OO:__setitem__", &index, &value)) return NULL;
  if(structure_ass_subscript(_self, index, value) < 0) return NULL;
  Py_RETURN_NONE;
}
#   undef LADA_STARTTRY
#   undef LADA_ENDTRY
