#include "PyladaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL pylada_ce_ARRAY_API
#include <python/include_numpy.h>

#include <algorithm>

#include <errors/exceptions.h>
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
# define PyMODINIT_FUNC void
#endif

#include "productilj.h"

namespace Pylada
{
  namespace ce
  {
    template<class T> void _outer_sum_impl( Py_ssize_t const _n,
                                            PyArrayObject* _result,
                                            Py_ssize_t const *_radices, 
                                            PyArrayObject * const * const _first,
                                            PyArrayObject * const * const _last )
    {
      T *out = (T*) PyArray_DATA(_result);
      for(Py_ssize_t i(0); i < _n; ++i, ++out)
      {
        T dummy(1);
        Py_ssize_t remainder = i;
        Py_ssize_t j(0);
        Py_ssize_t const *radix(_radices);
        for(PyArrayObject *const *vec(_first); vec != _last; ++j, ++vec, ++radix)
        {
          Py_ssize_t const component = remainder / (*radix);
          remainder %= (*radix);
          dummy *= *(T*) PyArray_GETPTR1(*vec, component);
        }
        *out += dummy;
      }
    }

    //! Checks whether a numpy array is an integer. 
    PyObject* outer_sum(PyObject *_module, PyObject *_args)
    {
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N == 1) Py_RETURN_NONE;
#     ifdef PYLADA_DEBUG
        if(N == 0)
        {
          PYLADA_PYERROR(TypeError, "sum_outer expects at least one arguments");
          return NULL;
        }
        PyArrayObject *result = (PyArrayObject*)PyTuple_GET_ITEM(_args, 0);
        if(not PyArray_Check(result))
        {
          PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                        "numpy arrays. The first one isn't.");
          return NULL;
        }
        if(not PyArray_ISCONTIGUOUS(result))
        {
          PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                        "c-contiguous numpy arrays. The first one isn't.");
          return NULL;
        }
        if(not PyArray_ISWRITEABLE(result))
        {
          PYLADA_PYERROR( TypeError, "sum_outer expects its first argument to be "
                        "writeable.");
          return NULL;
        }
        if(not (PyArray_ISFLOAT(result) or PyArray_ISINTEGER(result)))
        {
          PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                        "float or integer arrays." );
          return NULL;
        }
        if(PyArray_NDIM(result) != N-1)
        {
          PYLADA_PYERROR( TypeError, "incorrect shape for first argument to "
                        "outer_sum." );
          return NULL;
        }
        int const type = ((PyArrayObject*)result)->descr->type_num;
        for(Py_ssize_t i(1); i < N; ++i)
        {
          PyArrayObject* o = (PyArrayObject*) PyTuple_GET_ITEM(_args, i);
          if(not PyArray_Check(o))
          {
            PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                          "numpy arrays.");
            return NULL;
          }
          if(not PyArray_ISCONTIGUOUS(result))
          {
            PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                          "c-contiguous numpy arrays.");
            return NULL;
          }
          if(o->descr->type_num != type)
          {
            PYLADA_PYERROR( TypeError, "sum_outer expects all its arguments to be "
                          "arrays of the same type (float, int..).");
            return NULL;
          }
          if(PyArray_DIM(result, i-1) != PyArray_DIM(o, 0))
          {
            PYLADA_PYERROR( TypeError, "Dimensions of first argument and following "
                          "vectors do not match in outer_sum." );
            return NULL;
          }
        }
#     else
        PyArrayObject *result = (PyArrayObject*)PyTuple_GET_ITEM(_args, 0);
        int const type = ((PyArrayObject*)result)->descr->type_num;
#     endif

      Py_ssize_t const ndim = PyArray_NDIM(result);
      Py_ssize_t n = 1;
      for(Py_ssize_t i(0); i < ndim; ++i)
        n *= PyArray_DIM(result, i);

      PyArrayObject **vecs = new PyArrayObject*[ndim];
      if(not vecs)
      {
         PYLADA_PYERROR(internal, "Could not allocate memory.");
         return NULL;
      }
      Py_ssize_t *radices = new Py_ssize_t[ndim];
      if(not vecs)
      {
         PYLADA_PYERROR(internal, "Could not allocate memory.");
         delete[] vecs;
         return NULL;
      }
      Py_ssize_t const tobyte = PyArray_STRIDE(result, ndim-1);
      for(Py_ssize_t i(0); i < ndim; ++i)
      {
        vecs[i] = (PyArrayObject*) PyTuple_GET_ITEM(_args, i+1);
        radices[i] = PyArray_STRIDE(result, i) / tobyte;
      }
#     ifdef PYLADA_IFTYPE
#       error PYLADA_IFTYPE already defined
#     endif
#     define PYLADA_IFTYPE(TYPENUM, TYPE)                                        \
        case TYPENUM:                                                          \
           _outer_sum_impl<TYPE>(n, result, radices, vecs, vecs+ndim);         \
           break;
      switch(type)
      {
        PYLADA_IFTYPE(NPY_FLOAT,      npy_float)
        PYLADA_IFTYPE(NPY_DOUBLE,     npy_double)
        PYLADA_IFTYPE(NPY_LONGDOUBLE, npy_longdouble)
        PYLADA_IFTYPE(NPY_INT,        npy_int)
        PYLADA_IFTYPE(NPY_UINT,       npy_uint)
        PYLADA_IFTYPE(NPY_LONG,       npy_long)
        PYLADA_IFTYPE(NPY_ULONG,      npy_ulong)
        PYLADA_IFTYPE(NPY_LONGLONG,   npy_longlong)
        PYLADA_IFTYPE(NPY_ULONGLONG,  npy_ulonglong)
        PYLADA_IFTYPE(NPY_BYTE,       npy_byte)
        PYLADA_IFTYPE(NPY_UBYTE,      npy_ubyte)
        PYLADA_IFTYPE(NPY_SHORT,      npy_short)
        PYLADA_IFTYPE(NPY_USHORT,     npy_ushort)
      }
#     undef PYLADA_WITH_DATA_TYPE
      delete[] vecs;
      delete[] radices;
      Py_RETURN_NONE;
    }

    //! Methods table for ce module.
    static PyMethodDef methods_table[] = {
        {"outer_sum",  outer_sum, METH_VARARGS,
         "Sums to input n-tensor the n-tensor created from n vectors.\n\n"
         "From the vectors :py:math:`s_0, s_1, ..., s_N`, this function creates\n"
         "the n-dimensional tensor\n"
         ":py:math:`[T^{(N)}]_{i_0, ..., i_N}=\\prod_js_j^{(i_j)}``.\n\n"
         "The first argument should be the output n-tensor. The other\n"
         "arguments should be the n vectors of the correct dimensionality. The\n"
         "shape of the first arguments should be, loosely, the concatenated\n" 
         "shapes of the vectors. If the vectors are [0, 0], [1, 1], [2, 2, 2],\n"
         "then the shape of the first argument should be (2, 2, 3).\n"
         "It is possible, if inefficient, to use outer_sum with a single vector.\n"
         "It merely sums the vector to the first argument.\n"
         "It is also possible to use outer_sum without any vectors, in which case\n"
         "it does nothing.\n"
         ".. warning::\n\n"
         "   Don't use this if you do not understand the source code.\n"
         "   Or at the very least, first try things out with Pylada compiled\n"
         "   in debug mode. Otherwise, the input arrays are note checked\n"
         "   for correctness. If used incorrectly, expect segfaults.\n" },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
  }
}

PyMODINIT_FUNC initcppwrappers(void) 
{
  char const doc[] =  "Wrapper around C++ cluster expansion methods.";
  PyObject* module = Py_InitModule3("cppwrappers", Pylada::ce::methods_table, doc);
  if(not module) return;
  import_array(); // needed for NumPy 

  if (PyType_Ready(Pylada::ce::productiljiterator_type()) < 0) return;
  Py_INCREF(Pylada::ce::productiljiterator_type());

  PyModule_AddObject(module, "ProductILJ", (PyObject *)Pylada::ce::productiljiterator_type());
}
