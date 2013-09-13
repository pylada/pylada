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

#ifndef PYLADA_VFF_NODE_DATA_H
#define PYLADA_VFF_NODE_DATA_H

#include "PyladaConfig.h"

#include <vector>

//! \def check_structure(object)
//!      Returns true if an object is a node or subtype.
#define PyNodeData_Check(object) PyObject_TypeCheck(object, Pylada::vff::node_type())
//! \def PyNodeData_CheckExact(object)
//!      Returns true if an object is a node.
#define PyNodeData_CheckExact(object) object->ob_type == Pylada::vff::node_type()
      

namespace Pylada 
{
  namespace vff
  {
    extern "C" 
    {
      class EdgeData;
      //! Describes a node in a first neighbor net.
      struct NodeData
      {
        PyObject_HEAD 
        //! Holds list of weak pointers.
        PyObject *weakreflist;
        //! Holds possible gradient object.
        PyObject *gradient;
        //! Holds reference to other bonds.
        std::vector<EdgeData*> bonds;
        //! Holds reference to other an atom.
        crystal::Atom center;
        //! Index of the atom in the structure.
        long index;
      };
      //! Creates a new node.
      NodeData* PyNodeData_New();
      // Creates a new structure with a given type.
      NodeData* PyNode_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs);
      //! Adds an edge between two bonds. 
      bool PyNode_AddEdge(NodeData* _a, NodeData* _b, math::rVector3d &_trans);
      // Returns pointer to node type.
      PyTypeObject* node_type();

      //! Returns pointer to bond iterator type.
      PyTypeObject* bonditerator_type();
      //! Returns pointer to single-counting bond iterator type.
      PyTypeObject* dcbonditerator_type();
      //! Returns pointer to angle iterator type.
      PyTypeObject* angleiterator_type();
    }
  } // namespace vff

} // namespace Pylada

#endif
