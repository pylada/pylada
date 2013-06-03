#include "PyladaConfig.h"
#include "FCMangle.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL pylada_ewald_ARRAY_API
#define NO_IMPORT_ARRAY
#include <python/include_numpy.h>

#define PYLADA_NO_IMPORT
#include <python/python.h>

#include "ewald.h"

extern "C" void FC_GLOBAL( ewaldf, EWALDF )
                (
                  const int *const,    // verbosity
                  double * const,      // Energy
                  double * const,      // forces (reduced)
                  double * const,      // forces (cartesian)
                  const double *const, // stress
                  const int *const,    // number of atoms
                  const double *const, // reduced atomic coordinates.
                  const double *const, // atomic charges
                  const double *const, // real space cutoff
                  const double *const, // cell vectors
                  const int *const,    // dimension of arrays.
                  int * const          // ERROR
                );
namespace Pylada
{
  namespace pcm
  {
    PyObject* ewald(PyObject *_module, PyObject* _args, PyObject* _kwargs)  
    {
      static char *kwlist[] = { const_cast<char*>("cell"),
                                const_cast<char*>("positions"),
                                const_cast<char*>("charges"),
                                const_cast<char*>("cutoff"), NULL };
      PyObject *cell = NULL;
      PyObject *positions = NULL;
      PyObject *charges = NULL;
      double cutoff = 0;
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "OOOd:ewals",
                                          kwlist, &cell, &positions,
                                          &charges, &cutoff) )
        return NULL;

      const int verbosity(0);
      double energy(0);
      int error;
      npy_intp dims[3] = {PyArray_DIM(positions, 0), 3, 6};
      int const nptype = python::numpy::type<double>::value;
      python::Object forces = PyArray_ZEROS(2, dims, nptype, 0);
      if(not forces) return NULL;
      python::Object cforces = PyArray_ZEROS(2, dims, nptype, 0);
      if(not cforces) return NULL;
      python::Object stress = PyArray_ZEROS(1, &dims[2], nptype, 0);
      if(not stress) return NULL;
      int const n0 = dims[0];

      FC_GLOBAL( ewaldf, EWALDF )
      (
        &verbosity,                                 // verbosity
        &energy,                                    // Energy
        (double*) PyArray_DATA(forces.borrowed()),  // forces (reduced)
        (double*) PyArray_DATA(cforces.borrowed()), // forces (cartesian)
        (double*) PyArray_DATA(stress.borrowed()),  // stress
        &n0,                                        // number of atoms
        (double*) PyArray_DATA(positions),          // reduced atomic coordinates.
        (double*) PyArray_DATA(charges),            // atomic charges
        &cutoff,                                    // g-space cutoff in Ry.
        (double*) PyArray_DATA(cell),               // cell vectors
        &n0,                                        // dimension of arrays.
        &error
      );
      if(error == 1)
      {
        PYLADA_PYERROR(internal, "Could not find optimal alpha for ewald summation.");
        return NULL;
      }
      python::Object result = PyTuple_New(4);
      if(not result) return NULL;
      PyObject *pyenergy = PyFloat_FromDouble(energy);
      if(not pyenergy) return NULL;
      PyTuple_SET_ITEM(result.borrowed(), 0, pyenergy);
      PyTuple_SET_ITEM(result.borrowed(), 1, forces.release());
      PyTuple_SET_ITEM(result.borrowed(), 2, cforces.release());
      PyTuple_SET_ITEM(result.borrowed(), 3, stress.release());
      return result.release();
    }
  } // namespace pcm
} // namespace Pylada
