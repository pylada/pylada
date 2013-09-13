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

#include <iostream>
#include <dlfcn.h>

#include <boost/python/def.hpp>

#include <crystal/structure.h>
#include <opt/mpi.h>


#include "escan.hpp"

#ifdef PYLADA_DO_ESCAN
//! \cond
extern "C"
{
  void FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
  void FC_GLOBAL_(getvlarg, GETVLARG)();
  void FC_GLOBAL_(iaga_just_call_escan, IAGA_JUST_CALL_ESCAN)();
}
//! \endcond

void set_mpi(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_GLOBAL_(iaga_set_mpi, IAGA_SET_MPI)(&__commF);
}

void just_call_escan(boost::mpi::communicator const &_c)
{
  Py_BEGIN_ALLOW_THREADS;
  set_mpi(_c);
  FC_GLOBAL_(iaga_just_call_escan, IAGA_JUST_CALL_ESCAN)();
  Py_END_ALLOW_THREADS;
}

void just_call_genpot(boost::mpi::communicator const &_c)
{
  Py_BEGIN_ALLOW_THREADS;
  set_mpi(_c);
  FC_GLOBAL_(getvlarg, GETVLARG)();
  Py_END_ALLOW_THREADS;
}
#else 
  void just_call_genpot(boost::python::object const &_c)
  {
    PyErr_SetString(PyExc_ImportError, "Escan not found during compilation.");
    boost::python::throw_error_already_set();
  }
  void just_call_escan(boost::python::object const &_c)
  {
    PyErr_SetString(PyExc_ImportError, "Escan not found during compilation.");
    boost::python::throw_error_already_set();
  }
#endif 

namespace Pylada
{
  namespace python
  {
    namespace bp = boost::python;
    void expose_escan()
    {
      bp::def( "_call_escan", &just_call_escan,
               bp::arg("comm"), 
               "Calls escan, accepts a boost.mpi.communicator." );
      bp::def( "_call_genpot", &just_call_genpot,
               bp::arg("comm"), 
               "Calls genpot, accepts a boost.mpi.communicator." );
    }
  } // namespace python
} // namespace Pylada
