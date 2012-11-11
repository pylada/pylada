namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns cell as a numpy array. 
      static PyObject* structure_getcell(StructureData *_self, void *closure)
         { return python::wrap_to_numpy(_self->cell, (PyObject*)_self); }
      //! Sets cell from a sequence of 3x3 numbers.
      static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure);
      // Returns the scale.
      static PyObject* structure_getscale(StructureData *_self, void *closure);
      //! Sets the scale from a number.
      static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure);
      //! Gets the volume of the structure
      static PyObject* structure_getvolume(StructureData *_self, void *_closure)
      {
        types::t_real const scale = math::PyQuantity_AsReal(_self->scale);
        types::t_real const result = std::abs(_self->cell.determinant() * std::pow(scale, 3));
        return math::PyQuantity_FromCWithTemplate(result, _self->scale);
      }
    }
  
    // Sets cell from a sequence of three numbers.
    static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete cell attribute.");
        return -1;
      }
      return python::convert_to_matrix(_value, _self->cell) ? 0: -1;
    }
  
    // Returns the scale of the structure.
    static PyObject* structure_getscale(StructureData *_self, void *closure)
    {
      Py_INCREF(_self->scale); 
      return _self->scale;
    }
    // Sets the scale of the structure from a number.
    static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL) 
      {
        LADA_PYERROR(TypeError, "Cannot delete scale attribute.");
        return -1;
      }
      PyObject *result = math::PyQuantity_FromPy(_value, _self->scale);
      if(not result) return -1;
      PyObject *dummy = _self->scale;
      _self->scale = result;
      Py_DECREF(dummy);
      return 0;
    }
  }
} 
