namespace LaDa
{
  namespace crystal
  {
    //! Returns transformation matrix as numpy array.
    static PyObject* hftransform_gettransform(PyHFTObject *_self, void *closure)
       { return python::wrap_to_numpy(_self->transform, (PyObject*)_self); }
    //! Returns periodicity quotient as numpy array.
    static PyObject* hftransform_getquotient(PyHFTObject *_self, void *closure)
       { return python::wrap_to_numpy(_self->quotient, (PyObject*)_self); }
    //! Number of unit-cells in the supercell.
    static PyObject* hftransform_getsize(PyHFTObject *_self, void *closure)
       { return PyLong_FromLong(   _self->quotient[0] 
                                 * _self->quotient[1]
                                 * _self->quotient[2] ); }
    
  }
} 
