// Special methods for zinc-blende VFF.
namespace zincblende
{
  namespace 
  {
    PyArrayObject *get_bond_params(PyObject *_self, PyObject *_a, PyObject *_b)
    {
      PyObject *name = PyString_FromString("__getitem__");
      if(not name) return NULL;
      PyObject *args = PyTuple_Pack(2, _a, _b);
      if(not args) { Py_DECREF(name); return NULL; }
      PyArrayObject *result = (PyArrayObject*)PyObject_CallMethodObjArgs(_self, name, args, NULL);
      Py_DECREF(name);
      Py_DECREF(args);
      if(not result) return NULL;
#     ifdef PYLADA_DEBUG
        if(not PyArray_Check(result)) 
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected bond-parameters in vff.");
          return NULL;
        }
        if(PyArray_TYPE(result) != python::numpy::type<types::t_real>::value)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected array kind for bond-parameters.");
          return NULL;
        }
        if(PyArray_NDIM(result) != 1)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected number of bond-parameters in vff.");
          return NULL;
        }
        if(PyArray_DIM(result, 0) != 6)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected number of bond-parameters in vff.");
          return NULL;
        }
#     endif
      return (PyArrayObject*) result;
    }
    PyArrayObject *get_angle_params(PyObject *_self, PyObject *_a, PyObject *_b, PyObject *_c)
    {
      PyObject *name = PyString_FromString("__getitem__");
      if(not name) return NULL;
      PyObject *args = PyTuple_Pack(3, _a, _b, _c);
      if(not args) { Py_DECREF(name); return NULL; }
      PyArrayObject *result = (PyArrayObject*)PyObject_CallMethodObjArgs(_self, name, args, NULL);
      Py_DECREF(name); 
      Py_DECREF(args); 
      if(not result) return NULL;
#     ifdef PYLADA_DEBUG
        if(not PyArray_Check(result)) 
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected bond-parameters in vff.");
          return NULL;
        }
        if(PyArray_TYPE(result) != python::numpy::type<types::t_real>::value)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected array kind for bond-parameters.");
          return NULL;
        }
        if(PyArray_NDIM(result) != 1)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected number of bond-parameters in vff.");
          return NULL;
        }
        if(PyArray_DIM(result, 0) != 7)
        {
          Py_DECREF(result); 
          PYLADA_PYERROR(TypeError, "Unexpected number of bond-parameters in vff.");
          return NULL;
        }
#     endif
      return (PyArrayObject*) result;
    }
    const types::t_real fac6 = std::pow(std::sqrt(3e0)/2e0, 6) / types::t_real(6*5*4*3*2*1); 
    const types::t_real fac5 = std::pow(std::sqrt(3e0)/2e0, 5) / types::t_real(5*4*3*2*1); 
    const types::t_real fac4 = std::pow(std::sqrt(3e0)/2e0, 4) / types::t_real(4*3*2*1); 
    const types::t_real fac3 = std::pow(std::sqrt(3e0)/2e0, 3) / types::t_real(3*2*1); 
    const types::t_real fac2 = std::pow(std::sqrt(3e0)/2e0, 2) / types::t_real(2*1); 
    const types::t_real gfac6 = std::pow(std::sqrt(3e0)/2e0, 6) / types::t_real(5*4*3*2*1); 
    const types::t_real gfac5 = std::pow(std::sqrt(3e0)/2e0, 5) / types::t_real(4*3*2*1); 
    const types::t_real gfac4 = std::pow(std::sqrt(3e0)/2e0, 4) / types::t_real(3*2*1); 
    const types::t_real gfac3 = std::pow(std::sqrt(3e0)/2e0, 3);
    const types::t_real gfac2 = std::pow(std::sqrt(3e0)/2e0, 2);
#   ifdef PYLADA_BONDLENGTH
#     error PYLADA_BONDLENGTH already defined
#   endif 
#   ifdef PYLADA_ALPHA2
#     error PYLADA_ALPHA2 already defined
#   endif 
#   ifdef PYLADA_ALPHA3
#     error PYLADA_ALPHA3 already defined
#   endif 
#   ifdef PYLADA_ALPHA4
#     error PYLADA_ALPHA4 already defined
#   endif 
#   ifdef PYLADA_ALPHA5
#     error PYLADA_ALPHA5 already defined
#   endif 
#   ifdef PYLADA_ALPHA6
#     error PYLADA_ALPHA6 already defined
#   endif 
#   define PYLADA_BONDLENGTH(ARRAY) *(types::t_real*)PyArray_GETPTR1(ARRAY, 0)
#   define PYLADA_ALPHA2(ARRAY) fac2 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 1))
#   define PYLADA_ALPHA3(ARRAY) fac3 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 2))
#   define PYLADA_ALPHA4(ARRAY) fac4 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 3))
#   define PYLADA_ALPHA5(ARRAY) fac5 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 4))
#   define PYLADA_ALPHA6(ARRAY) fac6 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 5))
#   ifdef PYLADA_GAMMA
#     error PYLADA_GAMMA already defined
#   endif 
#   ifdef PYLADA_SIGMA
#     error PYLADA_SIGMA already defined
#   endif 
#   ifdef PYLADA_BETA2
#     error PYLADA_BETA2 already defined
#   endif 
#   ifdef PYLADA_BETA3
#     error PYLADA_BETA3 already defined
#   endif 
#   ifdef PYLADA_BETA4
#     error PYLADA_BETA4 already defined
#   endif 
#   ifdef PYLADA_BETA5
#     error PYLADA_BETA5 already defined
#   endif 
#   ifdef PYLADA_BETA6
#     error PYLADA_BETA6 already defined
#   endif 
#   define PYLADA_GAMMA(ARRAY) *(types::t_real*)PyArray_GETPTR1(ARRAY, 0)
#   define PYLADA_SIGMA(ARRAY) *(types::t_real*)PyArray_GETPTR1(ARRAY, 1)
#   define PYLADA_BETA2(ARRAY) fac2 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 2))
#   define PYLADA_BETA3(ARRAY) fac3 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 3))
#   define PYLADA_BETA4(ARRAY) fac4 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 4))
#   define PYLADA_BETA5(ARRAY) fac5 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 5))
#   define PYLADA_BETA6(ARRAY) fac6 * (*(types::t_real*)PyArray_GETPTR1(ARRAY, 6))
    PyObject *energy(PyObject *_module, PyObject *_args)
    {
      PyObject *self, *tree, *pystructure;
      if(not PyArg_ParseTuple(_args, "OOO", &self, &pystructure, &tree))
        return NULL;
      if(not crystal::check_structure(pystructure))
      {
        PYLADA_PYERROR(TypeError, "vff-energy: second argument should be a structure.");
        return NULL;
      }
      crystal::PyStructureObject* const structure = (crystal::PyStructureObject*)pystructure;
      if(not PyList_Check(tree)) 
      {
        PYLADA_PYERROR(TypeError, "vff-energy: third argument should be a list of nodes.");
        return NULL;
      }
      types::t_real const scale = python::get_quantity(structure->scale, "angstrom");
      types::t_real const scale2 = scale * scale;
      types::t_real energy = 0e0;

      // first loop over each atom in the structure.
      Py_ssize_t const N = PySequence_Size(tree);
      for(Py_ssize_t i(0); i < N; ++i)
      {
        NodeData *pynode = (NodeData*)PyList_GET_ITEM(tree, i);
#       ifdef PYLADA_DEBUG
          if(not PyNodeData_Check(pynode))
          {
            PYLADA_PYERROR(TypeError, "vff-energy: exepected a list of nodes as input.");
            return NULL;
          }
#       endif

        // for each atom, loop over each bond.
        std::vector<EdgeData*> :: const_iterator i_bond = pynode->bonds.begin();
        std::vector<EdgeData*> :: const_iterator const i_bond_end = pynode->bonds.end();
        for(; i_bond != i_bond_end; ++i_bond)
        {
          EdgeData* const bond0 = *i_bond;

          // bond_direction: avoid double-counting bonds
          bool const bond_direction = bond0->a == pynode;
          math::rVector3d const vector0 = bond_direction ?
            math::rVector3d(bond0->b->center->pos - pynode->center->pos + structure->cell * bond0->translation):
            math::rVector3d(bond0->a->center->pos - pynode->center->pos - structure->cell * bond0->translation);
        
          // get bond parameters.
          PyArrayObject * const bond0_params = get_bond_params( self, 
                                                                bond0->a->center->type, 
                                                                bond0->b->center->type );
          if(not bond0_params) return NULL;
          
          // compute e0 
          types::t_real const bondlength0 = PYLADA_BONDLENGTH(bond0_params);
          types::t_real const lambda0 = vector0.dot(vector0) * scale2 / bondlength0 - bondlength0;
          if(bond_direction) 
            energy +=  lambda0 * lambda0 * ( PYLADA_ALPHA2(bond0_params)
                         + lambda0 * (PYLADA_ALPHA3(bond0_params) 
                           + lambda0 * (PYLADA_ALPHA4(bond0_params) 
                            + lambda0 * (PYLADA_ALPHA5(bond0_params) 
                             + lambda0 * (PYLADA_ALPHA6(bond0_params) )))));

          std::vector<EdgeData*> :: const_iterator i_angle = i_bond + 1;
          for(; i_angle != i_bond_end; ++i_angle)
          {
            EdgeData* const bond1 = *i_angle;
            bool const angle_direction = bond1->a == pynode;
            math::rVector3d const vector1 = angle_direction ?
              math::rVector3d(bond1->b->center->pos - pynode->center->pos + structure->cell * bond1->translation):
              math::rVector3d(bond1->a->center->pos - pynode->center->pos - structure->cell * bond1->translation);
            // get bond and angle parameters.
            PyArrayObject * const bond1_params = get_bond_params( self, 
                                                                  bond1->a->center->type, 
                                                                  bond1->b->center->type );
            if(not bond0_params) return NULL;
            PyArrayObject * const angle_params =
                get_angle_params( self,
                                  angle_direction? bond1->b->center->type: bond1->a->center->type,
                                  pynode->center->type,
                                  bond_direction? bond0->b->center->type: bond0->a->center->type );
            if(not angle_params) return NULL;
            types::t_real const bondlength1 = PYLADA_BONDLENGTH(bond1_params);
            types::t_real const lambda1 = vector1.dot(vector1) * scale2 / bondlength1 - bondlength1;
            types::t_real const mean_length = std::sqrt(bondlength0*bondlength1); 
            types::t_real const beta = vector0.dot(vector1) * scale2 / mean_length
                                       - mean_length * PYLADA_GAMMA(angle_params);

            // add angle energy.
            energy +=  beta * beta * ( PYLADA_BETA2(angle_params)
                          + beta * (PYLADA_BETA3(angle_params) 
                           + beta * (PYLADA_BETA4(angle_params) 
                            + beta * (PYLADA_BETA5(angle_params) 
                             + beta * (PYLADA_BETA6(angle_params) )))));
            // add bond-angle energy
            energy += beta * (lambda0 + lambda1) * fac2 * PYLADA_SIGMA(angle_params);
          }
        } // loop over bonds
      } // loop over nodes
      return PyFloat_FromDouble(energy);
    } // energy function
   


    PyObject *jacobian(PyObject *_module, PyObject *_args)
    {
      PyObject *self, *tree, *pystructure;
      if(not PyArg_ParseTuple(_args, "OOO", &self, &pystructure, &tree))
        return NULL;
      if(not crystal::check_structure(pystructure))
      {
        PYLADA_PYERROR(TypeError, "vff-energy: second argument should be a structure.");
        return NULL;
      }
      crystal::PyStructureObject* const structure = (crystal::PyStructureObject*)pystructure;
      if(not PyList_Check(tree)) 
      {
        PYLADA_PYERROR(TypeError, "vff-energy: third argument should be a list of nodes.");
        return NULL;
      }
      types::t_real const scale = python::get_quantity(structure->scale, "angstrom");
      types::t_real const scale2 = scale * scale;
      types::t_real energy = 0e0;
      math::rMatrix3d stress = math::rMatrix3d::Zero(); 
      npy_intp dimensions[2] = {structure->atoms.size(), 3};
      PyArrayObject *pyforces = (PyArrayObject*)
                                   PyArray_ZEROS( 2, dimensions, 
                                                  python::numpy::type<types::t_real>::value,
                                                  0 ); 
      if(not pyforces) return NULL;

      // first loop over each atom in the structure.
      Py_ssize_t const N = PySequence_Size(tree);
      for(Py_ssize_t i(0); i < N; ++i)
      {
        NodeData *pynode = (NodeData*)PyList_GET_ITEM(tree, i);
#       ifdef PYLADA_DEBUG
          if(not PyNodeData_Check(pynode))
          {
            PYLADA_PYERROR(TypeError, "vff-energy: exepected a list of nodes as input.");
            Py_DECREF(pyforces);
            return NULL;
          }
#       endif

        // for each atom, loop over each bond.
        std::vector<EdgeData*> :: const_iterator i_bond = pynode->bonds.begin();
        std::vector<EdgeData*> :: const_iterator const i_bond_end = pynode->bonds.end();
        for(; i_bond != i_bond_end; ++i_bond)
        {
          EdgeData* const bond0 = *i_bond;

          // bond_direction: avoid double-counting bonds
          bool const bond_direction = bond0->a == pynode;
          NodeData* const endpoint0 = bond_direction ? bond0->b: bond0->a;
          math::rVector3d const vector0 = bond_direction ?
            math::rVector3d(endpoint0->center->pos - pynode->center->pos + structure->cell * bond0->translation):
            math::rVector3d(endpoint0->center->pos - pynode->center->pos - structure->cell * bond0->translation);
        
          // get bond parameters.
          PyArrayObject * const bond0_params = get_bond_params( self, 
                                                                pynode->center->type, 
                                                                endpoint0->center->type );
          if(not bond0_params) {Py_DECREF(pyforces); return NULL;}
          
          // compute e0 
          types::t_real const bondlength0 = PYLADA_BONDLENGTH(bond0_params);
          types::t_real const lambda0 = vector0.dot(vector0) * scale2 / bondlength0 - bondlength0;
          types::t_real const e0grad0 = 2e0 * scale2 / bondlength0 * lambda0 * ( 
                                        2e0 * PYLADA_ALPHA2(bond0_params) 
                                         + lambda0 * (3e0 * PYLADA_ALPHA3(bond0_params) 
                                           + lambda0 * (4e0 * PYLADA_ALPHA4(bond0_params) 
                                            + lambda0 * (5e0 * PYLADA_ALPHA5(bond0_params) 
                                             + lambda0 * (6e0 * PYLADA_ALPHA6(bond0_params) )))));
          if(bond_direction)
          {
            math::rVector3d const hold = vector0 * e0grad0;
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 0) -= hold[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 1) -= hold[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 2) -= hold[2];
 
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 0) += hold[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 1) += hold[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 2) += hold[2];

            stress += (vector0 * vector0.transpose()) * e0grad0;
            energy +=  lambda0 * lambda0 * ( PYLADA_ALPHA2(bond0_params)
                         + lambda0 * (PYLADA_ALPHA3(bond0_params) 
                           + lambda0 * (PYLADA_ALPHA4(bond0_params) 
                            + lambda0 * (PYLADA_ALPHA5(bond0_params) 
                             + lambda0 * (PYLADA_ALPHA6(bond0_params) )))));
          }

          std::vector<EdgeData*> :: const_iterator i_angle = i_bond + 1;
          for(; i_angle != i_bond_end; ++i_angle)
          {
            EdgeData* const bond1 = *i_angle;
            bool const angle_direction = bond1->a == pynode;
            NodeData* const endpoint1 = angle_direction ? bond1->b: bond1->a;
            math::rVector3d const vector1 = angle_direction ?
              math::rVector3d(endpoint1->center->pos - pynode->center->pos + structure->cell * bond1->translation):
              math::rVector3d(endpoint1->center->pos - pynode->center->pos - structure->cell * bond1->translation);
            // get bond and angle parameters.
            PyArrayObject * const bond1_params = get_bond_params( self, 
                                                                  pynode->center->type, 
                                                                  endpoint1->center->type );
            if(not bond1_params) {Py_DECREF(pyforces); return NULL;}
            PyArrayObject * const angle_params =
                get_angle_params( self,
                                  endpoint0->center->type,
                                  pynode->center->type,
                                  endpoint1->center->type );
            if(not angle_params) {Py_DECREF(pyforces); return NULL;}
            types::t_real const bondlength1 = PYLADA_BONDLENGTH(bond1_params);
            types::t_real const lambda1 = vector1.dot(vector1) * scale2 / bondlength1 - bondlength1;
            types::t_real const mean_length = std::sqrt(bondlength0*bondlength1); 
            types::t_real const beta = vector0.dot(vector1) * scale2 / mean_length
                                       - mean_length * PYLADA_GAMMA(angle_params);
            types::t_real const e0grad1 = 2e0 * scale2 / bondlength1 * lambda1 * ( 
                                          2e0 * PYLADA_ALPHA2(bond1_params) 
                                           + lambda1 * (3e0 * PYLADA_ALPHA3(bond1_params) 
                                             + lambda1 * (4e0 * PYLADA_ALPHA4(bond1_params) 
                                              + lambda1 * (5e0 * PYLADA_ALPHA5(bond1_params) 
                                               + lambda1 * (6e0 * PYLADA_ALPHA6(bond1_params) )))));
            // add stress/forces from angle bending.
            types::t_real const e1grad = scale2 / mean_length * beta * ( 
                                         2e0 * PYLADA_BETA2(angle_params) 
                                          + beta * (3e0 * PYLADA_BETA3(angle_params) 
                                            + beta * (4e0 * PYLADA_BETA4(angle_params) 
                                             + beta * (5e0 * PYLADA_BETA5(angle_params) 
                                              + beta * (6e0 * PYLADA_BETA6(angle_params) )))));
            math::rVector3d const anglegrad0 = e1grad * vector0;
            math::rVector3d const anglegrad1 = e1grad * vector1;
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 0) -= anglegrad0[0] + anglegrad1[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 1) -= anglegrad0[1] + anglegrad1[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 2) -= anglegrad0[2] + anglegrad1[2];
 
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 0) += anglegrad1[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 1) += anglegrad1[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 2) += anglegrad1[2];
 
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 0) += anglegrad0[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 1) += anglegrad0[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 2) += anglegrad0[2];

            math::rMatrix3d const _nonsym = vector0 * vector1.transpose();
            math::rMatrix3d const nonsym = _nonsym + _nonsym.transpose();
            stress += nonsym * e1grad;

            // add angle energy.
            energy +=  beta * beta * ( PYLADA_BETA2(angle_params)
                          + beta * (PYLADA_BETA3(angle_params) 
                           + beta * (PYLADA_BETA4(angle_params) 
                            + beta * (PYLADA_BETA5(angle_params) 
                             + beta * (PYLADA_BETA6(angle_params) )))));
            // add bond-angle energy
            energy += beta * (lambda0 + lambda1) * fac2 * PYLADA_SIGMA(angle_params);

            // add forces.
            math::rVector3d const bagrad0 
                = ( (2e0 * beta / bondlength0) * vector0 
                    + ((lambda0 + lambda1) / mean_length) * vector1 ) 
                  * (scale2 * fac2 * PYLADA_SIGMA(angle_params));
            math::rVector3d const bagrad1 
                = ( (2e0 * beta / bondlength1) * vector1 
                    + ((lambda0 + lambda1) / mean_length) * vector0 ) 
                  * (scale2 * fac2 * PYLADA_SIGMA(angle_params));
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 0) -= bagrad0[0] + bagrad1[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 1) -= bagrad0[1] + bagrad1[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, pynode->index, 2) -= bagrad0[2] + bagrad1[2];
 
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 0) += bagrad0[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 1) += bagrad0[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint0->index, 2) += bagrad0[2];
 
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 0) += bagrad1[0];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 1) += bagrad1[1];
            *(types::t_real*) PyArray_GETPTR2(pyforces, endpoint1->index, 2) += bagrad1[2];

            // now add stress.
            stress += (scale2 * fac2 * PYLADA_SIGMA(angle_params)) 
                      * ( 2e0 * beta * ( vector0*vector0.transpose()/bondlength0 
                                         + vector1*vector1.transpose()/bondlength1 )
                          + (lambda0 + lambda1) / mean_length * nonsym ); 
          }
        } // loop over bonds
      } // loop over nodes
      PyObject *pyenergy = PyFloat_FromDouble(energy);
      if(not pyenergy) { Py_DECREF(pyforces); return NULL; }
      PyObject *pystress = python::numpy::wrap_to_numpy(stress);
      if(not pystress) { Py_DECREF(pyenergy); Py_DECREF(pyforces); return NULL; }
      PyObject *result = PyTuple_New(3);
      if(not result) { Py_DECREF(pyenergy); Py_DECREF(pyforces); Py_DECREF(pystress); return NULL; }
      PyTuple_SET_ITEM(result, 0, pyenergy);
      PyTuple_SET_ITEM(result, 1, pystress);
      PyTuple_SET_ITEM(result, 2, (PyObject*)pyforces);
      return result;
    } // energy function
#   undef PYLADA_BONDLENGTH
#   undef PYLADA_ALPHA2
#   undef PYLADA_ALPHA3
#   undef PYLADA_ALPHA4
#   undef PYLADA_ALPHA5
#   undef PYLADA_ALPHA6
#   undef PYLADA_GAMMA
#   undef PYLADA_SIGMA
#   undef PYLADA_BETA2
#   undef PYLADA_BETA3
#   undef PYLADA_BETA4
#   undef PYLADA_BETA5
#   undef PYLADA_BETA6
  } // anonymous namespace
} // zincblende

