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

#ifndef PYLADA_VFF_EDGE_DATA_H
#define PYLADA_VFF_EDGE_DATA_H

#include "PyladaConfig.h"

//! \def check_structure(object)
//!      Returns true if an object is an edge or subtype.
#define PyEdgeData_Check(object) PyObject_TypeCheck(object, Pylada::vff::edge_type())
//! \def PyNodeData_CheckExact(object)
//!      Returns true if an object is a edge.
#define PyEdgeData_CheckExact(object) object->ob_type == Pylada::vff::edge_type()
      

namespace Pylada 
{
  namespace vff
  {
    extern "C" 
    {
      struct NodeData;
      //! \brief Describes an edge between first-neighbor nodes.
      //! This object is not meant to be returned directly in python.
      //! Rather, it is exists as a python object only so that python can do
      //! the memory management.
      struct EdgeData
      {
        PyObject_HEAD 
        //! Whether a periodic translation is involved.
        bool do_translate;
        //! Translation.
        math::rVector3d translation;
        //! Node A of bond.
        NodeData* a; 
        //! Node B of bond.
        NodeData* b; 
      };
      //! Creates a new edge.
      EdgeData* PyEdge_New();
      //! Returns pointer to node type.
      PyTypeObject* edge_type();
    }
    //! Returns tuple with both a and b.
    PyObject* edge_to_tuple(EdgeData* _data);
    //! Returns tuple with only one of a and b.
    PyObject* edge_to_tuple(NodeData* _self, EdgeData* _data);

  } // namespace vff

} // namespace Pylada

#endif

