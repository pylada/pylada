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
#include "FCMangle.h"

#include <unistd.h>
#include <boost/lexical_cast.hpp>

#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_arg.hpp>

#include <python/numpy_types.h>
#include <math/eigen.h>

#include "wavefunctions.hpp"

#ifdef PYLADA_DO_ESCAN
# include <boost/mpi/communicator.hpp>

  extern "C"
  {
    void FC_GLOBAL_(escan_wfns_init, ESCAN_WFNS_INIT)( int*, char const*, double const *,
                                                       double const *, double const *,
                                                       double const *, int const *, MPI_Fint* );
    void FC_GLOBAL_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_DIMENSIONS)(int*, int*, int*);
    void FC_GLOBAL_(escan_wfns_read, ESCAN_WFNS_READ)
      (int*, int*, int*, int*, int const*, double*, double*, double*, int*);
    void FC_GLOBAL_(escan_wfns_read_nokram, ESCAN_WFNS_READ_NOKRAM)
      (int*, int*, int*, int*, int const*, double*, double*, double*);
    void FC_GLOBAL_(escan_get_nr, ESCAN_GET_NR)(int *);
    void FC_GLOBAL_(escan_get_mr_n, ESCAN_GET_MR_N)(int *);
    void FC_GLOBAL_(escan_get_n1_n2_n3, ESCAN_GET_N1_N2_N3)(int*, int*, int*);
    void FC_GLOBAL_(escan_get_cell, ESCAN_GET_CELL)(double *, double *, double *);
    void FC_GLOBAL_(d3fft_comp, D3FFT_COMP)(double*, double*, int*);
    void FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)();
  }

  namespace Pylada
  {
    namespace python
    {
      namespace bp = boost::python;
      namespace bm = boost::mpi;
   
      void escan_get_cell( math::rMatrix3d &_cell )
      {
        double a[9];
        FC_GLOBAL_(escan_get_cell, ESCAN_GET_CELL)(a, a+3, a+6);
        for(size_t i(0); i < 9; ++i) _cell(i%3, i/3) = a[i];
      }

      bp::object positions(bm::communicator const &_comm)
      {
        math::rMatrix3d mesh_;
        int max_, nr_, n1, n3_, n2_;
        escan_get_cell(mesh_);
        FC_GLOBAL_(escan_get_nr, ESCAN_GET_NR)(&nr_);
        FC_GLOBAL_(escan_get_n1_n2_n3, ESCAN_GET_N1_N2_N3)(&n1, &n2_, &n3_);
        mesh_.col(0) /= types::t_real(n1);
        mesh_.col(1) /= types::t_real(n2_);
        mesh_.col(2) /= types::t_real(n3_);
        max_ = nr_ / _comm.size();

        bp::object result = python::numpy::create_array<NPY_DOUBLE>(max_, 3, true);
        npy_intp *strides = python::numpy::get_strides_pointer(result);
        double * array_data = python::numpy::get_data_pointer<double*>(result);
        for(size_t i(0); i < max_; ++i)
        {
          int u(i + _comm.rank()*(nr_/_comm.size()));
          math::rVector3d const a(u/(n3_*n2_), (u%(n3_*n2_)) / n3_, (u%(n3_*n2_)) % n3_);
          math::rVector3d const b( mesh_ * a ); 
          for(size_t j(0); j < 3; ++j)
          {
            char *where = reinterpret_cast<char*>(array_data) + i*strides[0]+j*strides[1];
            *reinterpret_cast<double*>(where) = b(j);
          }
        }
        return result;
      }

      bp::tuple read_wavefunctions( bp::object const &_escan, 
                                    bp::object const &_indices, 
                                    bp::object const &_kpoint,
                                    double _scale, 
                                    bm::communicator const &_comm ) 
      {
        if (bp::len(_kpoint) != 3)
        {
          PyErr_SetString(PyExc_ValueError, "kpoint should be a sequence of three elements.");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        double const kpoint[3] = { bp::extract<double>(_kpoint[0]),
                                   bp::extract<double>(_kpoint[1]), 
                                   bp::extract<double>(_kpoint[2]) };
        int const N = bp::extract<int>( _escan.attr("nbstates") );
        double const smooth = bp::extract<double>( _escan.attr("smooth") );
        double const kinscal = bp::extract<double>( _escan.attr("kinetic_scaling") );
        bool const is_krammer = bp::extract<bool>( _escan.attr("is_krammer") );
        int const pottype = bp::extract<int>(_escan.attr("potential") + 1);
        std::vector<int> indices;
        // extract indices.
        if( bp::len(_indices) == 0 )
        {
          try
          {
            int j = bp::extract<int>(_indices);
            if( j < 0 ) j += N;
            if( j < 0 or j >= N)
            {
              PyErr_SetString(PyExc_IndexError, "Wavefunction index out-of-range.\n");
              bp::throw_error_already_set();
              return bp::tuple();
            }
            indices.push_back(j + 1); // fortran indices.
          }
          catch (...) 
          {
            PyErr_SetString(PyExc_ValueError, "Second argument should integer or sequence.\n");
            bp::throw_error_already_set();
            return bp::tuple();
          }
        }
        else 
        {
          for(size_t i(0); i < bp::len(_indices); ++i)
            try
            { 
              int j = bp::extract<int>(_indices[i]);
              if( j < 0 ) j += N;
              if( j < 0 or j >= N)
              {
                PyErr_SetString(PyExc_IndexError, "Wavefunction index out-of-range.\n");
                bp::throw_error_already_set();
                return bp::tuple();
              }
              indices.push_back(j + 1); // fortran indices.
            }
            catch (...) 
            {
              PyErr_SetString( PyExc_ValueError, 
                               "Second argument should integer or sequence of integers.\n" );
              bp::throw_error_already_set();
              return bp::tuple();
            }
        }

        // prepares to read wavefunctions
        std::string const orig = bp::extract<std::string>(_escan.attr("WAVECAR"))();
        int a(orig.size()), b(indices.size());
        MPI_Comm __commC = (MPI_Comm) ( _comm ) ;
        MPI_Fint __commF = MPI_Comm_c2f( __commC );
        FC_GLOBAL_(escan_wfns_init, ESCAN_WFNS_INIT)( &a, orig.c_str(), &_scale, 
                                                      &smooth, &kinscal, kpoint, &pottype, &__commF);
        // gets dimensions.
        int n0, n1(indices.size()), n2, g0, g1(3);
        FC_GLOBAL_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_dimensions)(&n0, &n2, &g0);
        
        // Creates numpy objects.
        bp::object wfns = python::numpy::create_array<NPY_CDOUBLE>(n0, n1, n2, true);
        bp::object gpoints = python::numpy::create_array<NPY_DOUBLE>(g0, g1, true);
        bp::object projs = python::numpy::create_array<NPY_DOUBLE>(g0, true);
        bp::object inverse; 
        if(is_krammer) inverse = python::numpy::create_array<NPY_INT>(g0, true);
        
        python::numpy::get_pyarray_pointer(wfns)->dimensions[0]
          = python::numpy::get_pyarray_pointer(gpoints)->dimensions[0];

        // finally reads wavefunctions
        if(is_krammer)
          FC_GLOBAL_(escan_wfns_read, ESCAN_WFNS_READ)
                  ( 
                    &n0, &n1, &n2, &g0,  // dimensions
                    &(indices[0]),      // indices_ to wavefunctions
                    // pointer to wavefunctions data
                    python::numpy::get_data_pointer<double *const>(wfns),    
                    // pointer to gpoint  data
                    python::numpy::get_data_pointer<double *const>(gpoints), 
                    // pointer to projector data (smooth cutoff)
                    python::numpy::get_data_pointer<double *const>(projs),   
                    // pointer to -G indices
                    python::numpy::get_data_pointer<int* const>(inverse)            
                  );
        else
          FC_GLOBAL_(escan_wfns_read_nokram, ESCAN_WFNS_READ_NOKRAM)
                  ( 
                    &n0, &n1, &n2, &g0,  // dimensions
                    &(indices[0]),      // indices_ to wavefunctions
                    // pointer to wavefunctions data
                    python::numpy::get_data_pointer<double *const>(wfns),    
                    // pointer to gpoint  data
                    python::numpy::get_data_pointer<double *const>(gpoints), 
                    // pointer to projector data (smooth cutoff)
                    python::numpy::get_data_pointer<double *const>(projs)
                  );
        // and cleanup fortran arrays.
        FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
        // returns a 4-tuple.
        return bp::make_tuple(wfns, gpoints, positions(_comm), projs, inverse);
      }



  //   bp::tuple to_realspace( bp::object const &_escan, bp::object const &_gwfns,
  //                           bm::communicator const &_comm )
  //   {
  //     boost::mpi::communicator world;
  //     MPI_Comm __commC = (MPI_Comm) ( _comm ) ;
  //     MPI_Fint __commF = MPI_Comm_c2f( __commC );
  //     std::string const orig = bp::extract<std::string>(_escan.attr("_INCAR"))()
  //                              + "."
  //                              + boost::lexical_cast<std::string>(world.rank());
  //     int a(orig.size());
  //     FC_GLOBAL_(escan_wfns_init, ESCAN_WFNS_INIT)(&a, orig.c_str(), &__commF);
  //     // sanity checks
  //     if(not python::numpy::check_is_complex_array(_gwfns)) return bp::tuple();
  //
  //     PyArrayObject *array = reinterpret_cast<PyArrayObject*>(_gwfns.ptr());
  //     if(array->strides[0] != 16) 
  //     {
  //       PyErr_SetString(PyExc_ValueError, "Argument is not a contiguous "
  //                                         "numpy array of complexes. \n");
  //       bp::throw_error_already_set();
  //       return bp::tuple();
  //     }
  //     int n0, n2, g0;
  //     FC_GLOBAL_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_DIMENSIONS)(&n0, &n2, &g0);
  //     if( array->nd == 1 )
  //     {
  //       if(array->dimensions[0] != n0)
  //       {
  //         PyErr_SetString(PyExc_ValueError, "Unexpected array size. \n");
  //         bp::throw_error_already_set();
  //         FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //         return bp::tuple();
  //       }
  //     }
  //     else if( array->strides[0] * n0 < array->strides[1] ) 
  //     {
  //       PyErr_SetString(PyExc_ValueError, "Unexpected array size along first dimension. \n");
  //       bp::throw_error_already_set();
  //       FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //       return bp::tuple();
  //     }
  //    
  //     // creates output array
  //     FC_GLOBAL_(escan_get_mr_n, ESCAN_GET_MR_N)(&n0);
  //     if(n0 % 2 != 0) 
  //     {
  //       PyErr_SetString(PyExc_RuntimeError, "Wrong array dimension in pescan.");
  //       bp::throw_error_already_set();
  //       FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //       return bp::tuple();
  //     }
  //     std::vector<npy_intp> dims(1, n0>>1);
  //     for(int i(1); i < array->nd; ++i) dims.push_back(array->dimensions[i]);
  //     PyObject *result = PyArray_ZEROS(dims.size(), &dims[0], NPY_CDOUBLE, 1);
  //     if( result == NULL or PyErr_Occurred() != NULL )
  //     {
  //       bp::throw_error_already_set();
  //       FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //       return bp::tuple();
  //     }
  //     bp::object resob = bp::object(bp::handle<>(result));
  //
  //     // now loops through arrays and perform fft
  //     if( array->nd == 1 ) 
  //     {
  //       int sign = -1;
  //       FC_GLOBAL_(d3fft_comp, D3FFT_COMP)
  //       (
  //         (double*)array->data, 
  //         (double*)reinterpret_cast<PyArrayObject*>(result)->data,
  //         &sign
  //       );
  //     } 
  //     else
  //     {
  //       int dim = 0;
  //       PyArrayIterObject *real_iter 
  //         = reinterpret_cast<PyArrayIterObject*>(PyArray_IterAllButAxis(result, &dim));
  //       if( not real_iter )
  //       {
  //         PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
  //         bp::throw_error_already_set();
  //         FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //         return bp::tuple();
  //       }
  //       PyArrayIterObject *recip_iter 
  //         = reinterpret_cast<PyArrayIterObject*>
  //           (
  //             PyArray_IterAllButAxis(reinterpret_cast<PyObject*>(array), &dim)
  //           );
  //       // object should release on destruction
  //       bp::object dummyA(bp::handle<>(reinterpret_cast<PyObject*>(real_iter)));
  //       if( not recip_iter )
  //       {
  //         PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
  //         bp::throw_error_already_set();
  //         FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //         return bp::tuple();
  //       }
  //       // object should release on destruction
  //       bp::object dummyB(bp::handle<>(reinterpret_cast<PyObject*>(recip_iter)));
  //       while (real_iter->index < real_iter->size and recip_iter->index < recip_iter->size)  
  //       {
  //         int sign = -1;
  //         FC_GLOBAL_(d3fft_comp, D3FFT_COMP)( (double*)recip_iter->dataptr,
  //                                           (double*)real_iter->dataptr, &sign);
  //         PyArray_ITER_NEXT(real_iter);
  //         PyArray_ITER_NEXT(recip_iter);
  //       }
  //     }
  //     
  //     int nr;
  //     FC_GLOBAL_(escan_get_nr, ESCAN_GET_NR)(&nr);
  //     reinterpret_cast<PyArrayObject*>(result)->dimensions[0] = nr / _comm.size();
  //     FC_GLOBAL_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); 
  //     // finally creates real space position data.
  //     return bp::make_tuple(resob, positions(_comm));
  //   }
   
      void expose_wfns()
      {
        import_array();
        bp::def
        (
          "read_wavefunctions", 
          &read_wavefunctions,
          (bp::arg("escan"), bp::arg("indices"), bp::arg("kpoint"), bp::arg("scale"), bp::arg("comm")),
          "Context with temporary arrays to wavefunctions and corresponding g-vectors.\n\n"
          ":Parameters:\n"
          "  escan\n    Escan functional with which calculation were performed.\n"
          "  indices\n    index or indices for which to recover the wavefunctions. "
            "Indices of wavefunctions correspond to the eigenvalues "
            "returned by the functional during calculation.\n"
          "  comm\n    Communicator of same size as for calculations.\n"
          "  iskrammer : bool\n    True if escan uses krammer defeneracy.\n\n"
          ":return: (wavefunctions, g-vectors, projs, inverse).\n"
            "  - an spin by N by x matrix holding the N wavefuntions/spinor.\n",
            "  - a 3 by x matrix with each row a G-vector.\n"
            "  - a 3 by x matrix with each row a R-vector.\n"
            "  - one-dimensional array of real coefficients to smooth higher energy G-vectors.\n"
            "  - one-dimensional array of integer indices to map G-vectors to -G.\n"
        );
  //     bp::def( "to_realspace", &to_realspace,
  //              (bp::arg("escan"), bp::arg("wfns"), bp::arg("comm")),
  //              "Returns wavefunctions in real-space.\n\n"
  //              ":Parameters:\n"
  //              "  escan\n    Escan functional.\n"
  //              "  wfns : numpy array\n    Array of reciprocal-space wavefunctions.\n"
  //              "  comm\n    boost mpi communicator.\n"
  //              ":return: (numpy array real-space wfns, generator/array of positions\n");
      }


    } // namespace python
  } // namespace Pylada
#else 
  namespace Pylada
  {
    namespace python
    {
      namespace bp = boost::python;

      bp::tuple read_wavefunctions( bp::object const &_escan, 
                                    bp::object const &_indices, 
                                    bp::object const &_comm,
                                    bool _is_krammer)
      {
        PyErr_SetString(PyExc_ImportError, "Escan not found during compilation.");
        bp::throw_error_already_set();
        return bp::make_tuple();
      }

      bp::tuple to_realspace( bp::object const &_escan, bp::object const &_gwfns,
                              bp::object const &_comm )
      {
        PyErr_SetString(PyExc_ImportError, "Escan not found during compilation.");
        bp::throw_error_already_set();
        return bp::make_tuple();
      }

      void expose_wfns()
      {
        import_array();
        bp::def
        (
          "read_wavefunctions", 
          &read_wavefunctions,
          (bp::arg("escan"), bp::arg("indices"), bp::arg("comm"), bp::arg("iskrammer")),
          "Context with temporary arrays to wavefunctions and corresponding g-vectors.\n\n"
          ":Parameters:\n"
          "  escan\n    Escan functional with which calculation were performed.\n"
          "  indices\n    index or indices for which to recover the wavefunctions. "
            "Indices of wavefunctions correspond to the eigenvalues "
            "returned by the functional during calculation.\n"
          "  comm\n    Communicator of same size as for calculations.\n"
          "  iskrammer : bool\n    True if escan uses krammer defeneracy.\n\n"
          ":return: (wavefunctions, g-vectors, projs, inverse).\n"
            "  - an spin by N by x matrix holding the N wavefuntions/spinor.\n",
            "  - a 3 by x matrix with each row a G-vector.\n"
            "  - a 3 by x matrix with each row a R-vector.\n"
            "  - one-dimensional array of real coefficients to smooth higher energy G-vectors.\n"
            "  - one-dimensional array of integer indices to map G-vectors to -G.\n"
        );
        bp::def( "to_realspace", &to_realspace,
                 (bp::arg("escan"), bp::arg("wfns"), bp::arg("comm")),
                 "Returns wavefunctions in real-space.\n\n"
                 ":Parameters:\n"
                 "  escan\n    Escan functional.\n"
                 "  wfns : numpy array\n    Array of reciprocal-space wavefunctions.\n"
                 "  comm\n    boost mpi communicator.\n"
                 ":return: (numpy array real-space wfns, generator/array of positions\n");
      }
    }
  }
#endif
