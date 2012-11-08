#include "LaDaConfig.h"

#include <Python.h>
#include <structmember.h>
#define PY_ARRAY_UNIQUE_SYMBOL lada_math_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include <python/exceptions.h>
#include <python/wrap_numpy.h>
#include <math/misc.h>

#include "pybase.h"
#include "../edge/pybase.h"
#include "iterator.hpp"

#include "sequence.hpp"

// include private stuff from crystal so we can access pos and type faster.
#include "../../crystal/atom/getset.hpp"

namespace LaDa
{
  namespace vff
  {
    // Creates a new node.
    NodeData* PyNode_New()
    {
      NodeData* result = (NodeData*) node_type()->tp_alloc(node_type(), 0);
      if(not result) return NULL;
      result->weakreflist = NULL;
      new(&result->bonds) std::vector<EdgeData*>;
      new(&result->center) crystal::Atom;
      return result;
    }
    // Creates a new structure with a given type.
    NodeData* PyNode_NewWithArgs(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
      { return PyNode_New(); }

    static int traverse(NodeData *self, visitproc visit, void *arg)
    {
      Py_VISIT(self->center.borrowed());
      std::vector<EdgeData*>::const_iterator i_first = self->bonds.begin();
      std::vector<EdgeData*>::const_iterator const i_end = self->bonds.end();
      for(; i_first != i_end; ++i_first) Py_VISIT(*i_first);
      return 0;
    }
  
    static int gcclear(NodeData *self)
    { 
      self->center.release();
      std::vector<EdgeData*>::iterator i_first = self->bonds.begin();
      std::vector<EdgeData*>::iterator const i_end = self->bonds.end();
      for(; i_first != i_end; ++i_first) Py_CLEAR(*i_first);
      return 0;
    }

    // Function to deallocate a string atom.
    static void dealloc(NodeData *_self)
    {
      gcclear(_self);
  
      // calls destructor explicitely.
      PyTypeObject* ob_type = _self->ob_type;
      _self->~NodeData();

      ob_type->tp_free((PyObject*)_self);
    }
  
    // \brief Adds an edge between two bonds. 
    // \details Adds a bond to both nodes if it does not already exist.
    // \returns false and sets a python error if something bad occured.
    bool PyNode_AddEdge(NodeData* _a, NodeData* _b, math::rVector3d &_trans)
    {
      if(_a == _b and math::is_null(_trans))
      {
        LADA_PYERROR(ValueError, "Cannot add link from one atom to itself (same periodic image).");
        return false;
      }
      // check if edge already exists.
      // We need only check edges from one node.
      std::vector<EdgeData*>::const_iterator i_first = _a->bonds.begin();
      std::vector<EdgeData*>::const_iterator const i_end = _a->bonds.end();
      for(; i_first != i_end; ++i_first)
        if((*i_first)->a == _a and (*i_first)->b == _b)
        {
          if(math::eq((*i_first)->translation, _trans)) return true;
        }
        else if((*i_first)->a == _b and (*i_first)->b == _a) 
        {
          if(math::eq((*i_first)->translation, -_trans)) return true;
        }
      EdgeData *bond =  PyEdge_New();
      if(not bond) return false;
      bond->a = _a;
      Py_INCREF(_a);
      bond->b = _b;
      if(bond->a != bond->b) Py_INCREF(_b);
      bond->translation = _trans;
      bond->do_translate = not math::is_null(_trans);
      _a->bonds.push_back(bond);
      if(_a != _b)
      { 
        _b->bonds.push_back(bond);
        Py_INCREF(bond);
      }
      return true;
    }

    static PyObject* add_edge(NodeData* _self, PyObject* _args, PyObject *_kwargs)
    {
      PyObject *endpoint = NULL;
      PyObject *_trans = NULL;
      static char *kwlist[] = { const_cast<char*>("endpoint"),
                                const_cast<char*>("translation"), NULL };
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "O|O:add_edge",
                                          kwlist, &endpoint, &_trans) )
        return NULL;

      if(not PyNodeData_Check(endpoint)) 
      {
        LADA_PYERROR(TypeError, "add_edge expects a NodeData end-point.");
        return NULL;
      }
      math::rVector3d translation(0,0,0);
      if(_trans and not python::convert_to_vector(_trans, translation)) return NULL;
      
      if(not PyNode_AddEdge(_self, (NodeData*)endpoint, translation)) Py_RETURN_FALSE;
      Py_RETURN_TRUE;
    }

    static PyObject* getpos(NodeData* _self, void *closure)
      { return crystal::lada_atom_getpos((crystal::PyAtomObject*)_self->center.borrowed(), closure); }
    static int setpos(NodeData* _self, PyObject* _value, void *closure)
      { return crystal::lada_atom_setpos((crystal::PyAtomObject*)_self->center.borrowed(), _value, closure); }
    static PyObject* gettype(NodeData* _self, void *closure)
      { return crystal::lada_atom_gettype((crystal::PyAtomObject*)_self->center.borrowed(), closure); }
    static int settype(NodeData* _self, PyObject* _value, void *closure)
      { return crystal::lada_atom_settype((crystal::PyAtomObject*)_self->center.borrowed(), _value, closure); }
    static PyObject* getcenter(NodeData* _self, void *closure)
      { return _self->center.new_ref(); }

    static int init(PyObject* _self, PyObject *_args, PyObject *_kwargs)
    {
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N != 1)
      {
        LADA_PYERROR(TypeError, "Node expects an atom as its only and required argument." );
        return -1;
      }
      if(_kwargs and PyDict_Size(_kwargs) != 0)
      {
        LADA_PYERROR(TypeError, "Node does not take keyword arguments.");
        return -1;
      }
      PyObject* node = PyTuple_GET_ITEM(_args, 0);
      if(not PyAtom_Check(node))
      {
        LADA_PYERROR(TypeError, "Node expects an atom as its only and required argument." );
        return -1;
      }
      ((NodeData*)_self)->center = crystal::Atom::acquire(node);
      return 0;
    }

    static PyObject* make_sure_itsdefined(NodeData* _in) 
      { return bond_iterator_create<&bonditerator_type>(_in); }
    static PyObject* make_sure_itsdefined2(NodeData* _in) 
      { return bond_iterator_create<&dcbonditerator_type>(_in); }

    // Returns pointer to node type.
    PyTypeObject* node_type()
    {
      static PyMappingMethods as_mapping = {
          (lenfunc)size,
          (binaryfunc)subscript,
          (objobjargproc)NULL
      };
      static PySequenceMethods as_sequence = {
          (lenfunc)size,                              /* sq_length */
          (binaryfunc)NULL,                           /* sq_concat */
          (ssizeargfunc)NULL,                         /* sq_repeat */
          (ssizeargfunc)getitem,                      /* sq_item */
          (ssizessizeargfunc)NULL,                    /* sq_slice */
          (ssizeobjargproc)NULL,                      /* sq_ass_item */
          (ssizessizeobjargproc)NULL,                 /* sq_ass_slice */
          (objobjproc)NULL,                           /* sq_contains */
          (binaryfunc)NULL,                           /* sq_inplace_concat */
          (ssizeargfunc)NULL,                         /* sq_inplace_repeat */
      };
#     ifdef LADA_DECLARE
#       error LADA_DECLARE already defined.
#     endif
#     define LADA_DECLARE(name, doc) \
        { const_cast<char*>(#name), (getter) get ## name, \
          (setter) set ## name, const_cast<char*>(doc) }
      
      static PyGetSetDef getsetters[] = {
          LADA_DECLARE(pos,  "Alias to wrapped atom's position."),
          LADA_DECLARE(type, "Alias to wrapped atom's type."),
          { const_cast<char*>("center"), (getter) getcenter, 
            NULL, const_cast<char*>("Wrapped atom.") },
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
#     define LADA_DECLARE(name, func, args, doc) \
        {#name, (PyCFunction)func, METH_ ## args, doc} 
      static PyMethodDef methods[] = {
          LADA_DECLARE( link,  link, KEYWORDS,
                        "Adds a bond between to nodes."),
          LADA_DECLARE(clear,  clear, NOARGS, "Removes all bonds."),
          LADA_DECLARE(sc_bond_iter, 
                       bond_iterator_create<dcbonditerator_type>, 
                       NOARGS, 
                       "Iterates over bonds without double counting.\n\n"
                       "This function allows iterating over the bonds of a"
                       "structure without double counting. Nevertheless, it\n"
                       "proceed via a nested loop:\n\n"
                       ">>> for node in structure_tree:\n"
                       "...   for other, vector in node.sc_bond_iter():\n"
                       "...     perimag = dot(structure.cell, vector)\n"
                       "...     bond_vector = perimag + other.pos - node.pos\n\n"
                       "The above loops first over the atoms in a structure"
                       "tree obtained from :py:meth:`Functional.build_tree`,\n"
                       "and then over the bonds for that center."
                       "The resulting ``bond_vector`` is the vector between the\n"
                       "atom (in the first loop) to the end-point of the bond."),
          LADA_DECLARE( angle_iter, angle_iterator_create, NOARGS, 
                        "Iterates over angles centered on this node." ),
          {NULL}  /* Sentinel */
      };
#     undef LADA_DECLARE
      static PyTypeObject dummy = {
          PyObject_HEAD_INIT(NULL)
          0,                                 /*ob_size*/
          "lada.vff.cppwrappers.Node",   /*tp_name*/
          sizeof(NodeData),                  /*tp_basicsize*/
          0,                                 /*tp_itemsize*/
          (destructor)dealloc,               /*tp_dealloc*/
          0,                                 /*tp_print*/
          0,                                 /*tp_getattr*/
          0,                                 /*tp_setattr*/
          0,                                 /*tp_compare*/
          0,                                 /*tp_repr*/
          0,                                 /*tp_as_number*/
          &as_sequence,                      /*tp_as_sequence*/
          &as_mapping,                       /*tp_as_mapping*/
          0,                                 /*tp_hash */
          0,                                 /*tp_call*/
          0,                                 /*tp_str*/
          0,                                 /*tp_getattro*/
          0,                                 /*tp_setattro*/
          0,                                 /*tp_as_buffer*/
          Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_HAVE_ITER, /*tp_flags*/
          "Defines a node.\n\n"              /*tp_doc*/
          "Node of a first neighbor tree. It wraps around an atom and allows easy\n"
          "iteration on bonds and angles (once bonds have been added).\n",
          (traverseproc)traverse,            /* tp_traverse */
          (inquiry)gcclear,                  /* tp_clear */
          0,		                             /* tp_richcompare */
          0,                                 /* tp_weaklistoffset */
          (getiterfunc)&bond_iterator_create<bonditerator_type>,/* tp_iter */
          0,		                             /* tp_iternext */
          methods,                           /* tp_methods */
          0,                                 /* tp_members */
          getsetters,                        /* tp_getset */
          0,                                 /* tp_base */
          0,                                 /* tp_dict */
          0,                                 /* tp_descr_get */
          0,                                 /* tp_descr_set */
          0,                                 /* tp_dictoffset */
          init,                              /* tp_init */
          0,                                 /* tp_alloc */
          (newfunc)PyNode_NewWithArgs,       /* tp_new */
      };
      return &dummy;
    }

    // Returns bond iterator type, with double counting.
    PyTypeObject* bonditerator_type()
    { 
      static PyTypeObject type = {
          PyObject_HEAD_INIT(NULL)
          0,                                          /*ob_size*/
          "lada.vff.cppwrappers.BondIterator",        /*tp_name*/
          sizeof(BondIterator),                       /*tp_basicsize*/
          0,                                          /*tp_itemsize*/
          (destructor)bond_iterator_dealloc,          /*tp_dealloc*/
          0,                                          /*tp_print*/
          0,                                          /*tp_getattr*/
          0,                                          /*tp_setattr*/
          0,                                          /*tp_compare*/
          0,                                          /*tp_repr*/
          0,                                          
          0,                                          /*tp_as_sequence*/
          0,                                          /*tp_as_mapping*/
          0,                                          /*tp_hash */
          0,                                          /*tp_call*/
          0,                                          /*tp_str*/
          0,                                          /*tp_getattro*/
          0,                                          /*tp_setattro*/
          0,                                          /*tp_as_buffer*/
          Py_TPFLAGS_HAVE_ITER,
          "Iterator over bonds.\n\n"
          "Yields a tuple ``(center, v)``, where ``center``\n"
          "is the  nearest-neighbor :py:class:`Node` at the other\n"
          "end of the bond, and ``v`` is  a vector in fractional\n"
          "coordinates from which the vector linking the two end\n"
          "points can be obtained:\n\n"
          ">>> for center, v in node: \n"
          ">>>    atob = center.pos + dot(structure.cell, v) - node.pos\n\n"
          "The code above loops over all the bonds to a given node. The\n"
          "vector ``atob`` links the node to the correct periodic image\n"
          "of the node ``center`` at the end of the bond.",
          0,                                          /* tp_traverse */
          0,                                          /* tp_clear */
          0,                                          /* tp_richcompare */
          0,		                                      /* tp_weaklistoffset */
          (getiterfunc)get_self,                      /* tp_iter */
          (iternextfunc)bond_iterator_next,           /* tp_iternext */
      };
      return &type;
    }

    // Returns bond iterator type, without double counting.
    PyTypeObject* dcbonditerator_type()
    { 
      static PyTypeObject type = {
          PyObject_HEAD_INIT(NULL)
          0,                                          /*ob_size*/
          "lada.vff.cppwrappers.DcBondIterator",      /*tp_name*/
          sizeof(BondIterator),                       /*tp_basicsize*/
          0,                                          /*tp_itemsize*/
          (destructor)bond_iterator_dealloc,          /*tp_dealloc*/
          0,                                          /*tp_print*/
          0,                                          /*tp_getattr*/
          0,                                          /*tp_setattr*/
          0,                                          /*tp_compare*/
          0,                                          /*tp_repr*/
          0,                                          
          0,                                          /*tp_as_sequence*/
          0,                                          /*tp_as_mapping*/
          0,                                          /*tp_hash */
          0,                                          /*tp_call*/
          0,                                          /*tp_str*/
          0,                                          /*tp_getattro*/
          0,                                          /*tp_setattro*/
          0,                                          /*tp_as_buffer*/
          Py_TPFLAGS_HAVE_ITER,
          "Iterator over bonds, without double counting."
          "Yields a tuple ``(center, v)``, where ``center``\n"
          "is the  nearest-neighbor :py:class:`Node` at the other\n"
          "end of the bond, and ``v`` is  a vector in fractional\n"
          "coordinates from which the vector linking the two end\n"
          "points can be obtained:\n\n"
          ">>> from lada.vff import build_tree\n"
          ">>> tree = build_tree(structure)\n"
          ">>> for node in tree: \n"
          ">>>   for center, v in node.sc_bond_iter(): \n"
          ">>>      atob = center.pos + dot(structure.cell, v) - node.pos\n\n"
          "The code above loops over *all* the bonds in a structure, *without*\n"
          "double counting. The  vector ``atob`` links the node to the correct\n"
          "periodic image of the node ``center`` at the end of the bond.\n",
          0,                                          /* tp_traverse */
          0,                                          /* tp_clear */
          0,                                          /* tp_richcompare */
          0,		                              /* tp_weaklistoffset */
          (getiterfunc)get_self,                      /* tp_iter */
          (iternextfunc)dcbond_iterator_next,           /* tp_iternext */
      };
      return &type;
    }

    // Returns angle iterator type
    PyTypeObject* angleiterator_type()
    { 
      static PyTypeObject type = {
          PyObject_HEAD_INIT(NULL)
          0,                                          /*ob_size*/
          "lada.vff.cppwrappers.AngleIterator",       /*tp_name*/
          sizeof(AngleIterator),                      /*tp_basicsize*/
          0,                                          /*tp_itemsize*/
          (destructor)angle_iterator_dealloc,         /*tp_dealloc*/
          0,                                          /*tp_print*/
          0,                                          /*tp_getattr*/
          0,                                          /*tp_setattr*/
          0,                                          /*tp_compare*/
          0,                                          /*tp_repr*/
          0,                                          
          0,                                          /*tp_as_sequence*/
          0,                                          /*tp_as_mapping*/
          0,                                          /*tp_hash */
          0,                                          /*tp_call*/
          0,                                          /*tp_str*/
          0,                                          /*tp_getattro*/
          0,                                          /*tp_setattro*/
          0,                                          /*tp_as_buffer*/
          Py_TPFLAGS_HAVE_ITER,
          "Iterator over angles.",
          0,                                          /* tp_traverse */
          0,                                          /* tp_clear */
          0,                                          /* tp_richcompare */
          0,		                                      /* tp_weaklistoffset */
          (getiterfunc)get_self,                      /* tp_iter */
          (iternextfunc)angle_iterator_next,          /* tp_iternext */
      };
      return &type;
    }

  } // namespace vff
} // namespace LaDa
