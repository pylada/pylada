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

#ifndef PYLADA_ENUM_NDIMITERATOR_H
#define PYLADA_ENUM_NDIMITERATOR_H
#include "PyladaConfig.h"

#include <Python.h>

#include <vector>

//! \def PyNDimIterator_Check(object)
//!      Returns true if an object is a struture or subtype.
#define PyNDimIterator_Check(object) PyObject_TypeCheck(object, Pylada::crystal::structure_type())
//! \def PyNDimIterator_CheckExact(object)
//!      Returns true if an object is a structure.
#define PyNDimIterator_CheckExact(object) object->ob_type == Pylada::crystal::structure_type()
      

namespace Pylada
{
  namespace ce
  {
    extern "C" 
    {
      //! Product iterator over a sequence, with i < j < ...
      struct ProductILJIterator
      {
        PyObject_HEAD 
        //! Reference to sequence. 
        PyObject *sequence;
        //! Size of the sequence.
        Py_ssize_t N;
        //! Inner counter;
        std::vector<Py_ssize_t> counter;
        //! Whether this is the first iteration.
        bool is_first;
      };
      //! Creates a new structure.
      ProductILJIterator* PyNDimIterator_New();
      //! Creates a new structure with a given type.
      ProductILJIterator* PyNDimIterator_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Creates a new structure with a given type, also calling initialization.
      ProductILJIterator* PyNDimIterator_NewFromArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      // Returns pointer to structure type.
      PyTypeObject* productiljiterator_type();
    }
  }
}
#endif
