###############################
#  This file is part of PyLaDa.
#
#  Copyright (C) 2013 National Renewable Energy Lab
# 
#  PyLaDa is a high throughput computational platform for Physics. It aims to make it easier to submit
#  large numbers of jobs on supercomputers. It provides a python interface to physical input, such as
#  crystal structures, as well as to a number of DFT (VASP, CRYSTAL) and atomic potential programs. It
#  is able to organise and launch computational jobs on PBS and SLURM.
# 
#  PyLaDa is free software: you can redistribute it and/or modify it under the terms of the GNU General
#  Public License as published by the Free Software Foundation, either version 3 of the License, or (at
#  your option) any later version.
# 
#  PyLaDa is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
#  Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along with PyLaDa.  If not, see
#  <http://www.gnu.org/licenses/>.
###############################

class Vff(object): 
  """ A valence force field functional for zinc-blende. """
  def __init__(self): 
    super(Vff, self).__init__()

    self._parameters = {}
    """ Holds vff parameters. """

  def __getitem__(self, index):
    """ VFF parameter access. 
    
        VFF parameters are composed of bond and angle parameters.
        
        Bond parameters can be accessed as:

        >>> functional['Au', 'Pd']
        array([2.62332, 21.6739, -112.0, 150.0, 0, 0])

        The order of the two species does not matter:

        >>> functional['Au', 'Pd'] is functional['Pd', 'Au']
        True

        The first item of the return array is the bond-length. The rest are
        the bond stretching parameters (second to 6th order).

        Angle parameters can be accessed the same way:

        >>> functional['Au', 'Pd', 'Cu']
        (-0.3333, -4.099, 9.3703)

        The first argument is the cosine of the equilibrium angle.
        The order of the first and last specie does not matter. The center of
        the angle is the central (second) specie.
    """
    if isinstance(index, str): index = index.split('-')
    index = '-'.join(sorted(str(u) for u in index))
    return self._parameters[index]

  def __setitem__(self, index, value):
    """ Adds/sets VFF parameters.
    
        VFF parameters are composed of bond and angle parameters.
        
        Bond parameters can be added/set as:

        >>> functional['Au', 'Pd'] = 2.62332, 21.6739, -112.0, 150.0, 0, 0

        The order of the two species does not matter. At this point, units
        should be based on angstroms. This is not checked in any way.

        The first item of the return array is the bond-length. The rest are
        the bond stretching parameters (second to 6th order).

        Angle parameters can be added/set in the same way:

        >>> functional['Au', 'Pd', 'Cu'] = -0.3333, -4.099, 9.3703

        The first argument is the cosine of the equilibrium angle.
        In the case of a perfect tetrahedron, it can be set from the string
        "tet". Furthermore, it can be set from a quantity signed in degrees
        or radians (see quantity_).

        The rest of the parameters are the bond-stretching parameters from
        2nd to 6th order. Their units should be based on angstroms though
        this is not checked in anyway.
    """
    from numpy import cos, array
    from quantities import radian, meter, angstrom, newton
    from ..error import ValueError

    if isinstance(index, str): index = index.split('-')
    index = '-'.join(sorted(str(u) for u in index))
    # special case where the angle is given as "tet"
    value = list(value)
    maxsize = 7 if index.count('-') == 2 else 6
    if len(value) < 2 or len(value) > maxsize:
      raise ValueError( 'Expects no less than two and no more than 6 '         \
                        'parameters.')
    if maxsize == 7:
      if isinstance(value[0], str):
        if value[0][:3].lower() != 'tet':
          raise ValueError( 'If a string, the first argument to angle '        \
                            'parameters should be "tet". ')
        value[0] = -1e0/3e0
      elif hasattr(value[0], 'rescale'):
        value[0] = cos(value[0].rescale(radian).magnitude)
      if hasattr(value[1], 'rescale'):
        value[1] = float(value[1].rescale(newton/meter))
      for i in xrange(2, len(value)):
        if hasattr(value[i], 'rescale'):
          value[i] = float(value[i].rescale(newton/meter/angstrom**(i-2)))
    else:
      for i in xrange(1, len(value)):
        if hasattr(value[i], 'rescale'):
          value[i] = float(value[i].rescale(newton/meter/angstrom**(i-1)))
    value = array(value).flatten()
    
    self._parameters[index] = array( value.tolist()
                                     + [0]*(maxsize - len(value)) )

  def __contains__(self, index):
    """ True if bond or angle parameters are declared. """
    if isinstance(index, str): index = index.split('-')
    index = '-'.join(sorted(str(u) for u in index))
    return index in self._parameters

  def energy(self, structure, _tree=None):
    """ Evaluates energy alone. """
    from quantities import newton, angstrom, meter, eV
    from .cppwrappers import _zb_energy
    from . import build_tree

    # first, build tree
    tree = build_tree(structure) if _tree is None else _tree
    # then call cpp function.
    return _zb_energy(self, structure, tree)                                   \
           * (newton/meter*angstrom**2).rescale(eV)


  def jacobian(self, structure, _tree=None):
    """ Jacobian of the structure. 

        The jacobian is separated into a stress component and a atomic-gradient
        part. The stress is normalized by -1./volume, e.g. as per physical
        definition.
    """
    from numpy.linalg import det
    from quantities import angstrom, eV, meter, newton
    from .cppwrappers import _zb_jacobian
    from . import build_tree


    # first, build tree
    tree = build_tree(structure) if _tree is None else _tree
    energy, stress, forces = _zb_jacobian(self, structure, tree) 
    volume = det(structure.scale.rescale(angstrom)*structure.cell) * angstrom**3
    stress = -stress / volume * (newton / meter * angstrom * angstrom).rescale(eV)
    forces = forces * ( (newton / meter * angstrom).rescale(eV/angstrom)
                        / float(structure.scale.rescale(angstrom)) )
    return stress, forces

  def __call__(self, structure, _tree=None):
    """ Evaluates energy and forces on a structure. """
    from numpy.linalg import det
    from quantities import angstrom, eV, meter, newton
    from .cppwrappers import _zb_jacobian
    from . import build_tree


    # first, build tree
    tree = build_tree(structure) if _tree is None else _tree
    energy, stress, forces = _zb_jacobian(self, structure, tree) 
    volume = det(structure.scale.rescale(angstrom)*structure.cell) * angstrom**3

    result = structure.copy()
    result.energy = energy * (newton/meter*angstrom*angstrom).rescale(eV)
    result.stress = -stress / volume                                           \
                    * (newton / meter * angstrom * angstrom).rescale(eV)
    forceunits = (newton / meter * angstrom).rescale(eV/angstrom)              \
                 / float(structure.scale.rescale(angstrom)) 
    for atom, force in zip(result, forces):
      atom.gradient = force * forceunits
    return result

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly representation of the functional. """
    from ..tools.uirepr import template_ui_repr
    results = template_ui_repr(self, imports, name, defaults, ['_parameters'])
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    for key, value in self._parameters.iteritems():
      index = ', '.join(repr(u) for u in key.split('-'))
      key = "{0}[{1}]".format(name, index)
      results[key] = ', '.join(str(u) for u in value)
    return results

  def __repr__(self, defaults=True, name=None):
    """ Returns representation of this instance """
    from ..tools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)
    

  def _pyeval(self, structure):
    """ Evaluates energy and forces on a structure. 
    
        This follows a python implementation, rather than the faster c
        implementation.
    """
    from numpy import zeros
    from numpy.linalg import det
    from quantities import eV, angstrom
    from . import build_tree
    # creates result structure.
    result = structure.copy()
    for atom in result:
      atom.gradient = zeros(3, dtype='float64') * eV / angstrom
    result.stress = zeros((3,3), dtype='float64') * eV 
    result.energy = 0e0 * eV
   
    # creates tree and loop over structure.
    tree = build_tree(result)
    scale = float(result.scale.rescale(angstrom))
    for node in tree: 
      result.energy += self._evaluate_bonds( node, scale, result.cell,
                                             result.stress )
      result.energy += self._evaluate_angles( node, scale, result.cell,
                                              result.stress )
    result.stress *= -1e0/det(result.cell*result.scale) * angstrom**(-3)
    return result

  def _evaluate_bonds(self, node, scale, cell, stress=None):
    """ Evaluates bond-stretching modes.
    
        .. note:: stress is not normalized by -1/volume yet.
    """
    from numpy import dot, sum, array, sqrt, outer
    from quantities import newton, meter, angstrom, eV
    bondparams = array([1e0, 5e-1/sqrt(3e0), 1e0/16e0, sqrt(3e0)/8e0, 1e0/640])
    gbondparams = array([ 1.5e0, 3e0*sqrt(3e0)/8.0, 3e0/16e0,
                          3e0*sqrt(3e0)/128e0, 0.00703125 ])
    scale2 = scale * scale

    energy = 0
    stressunits = (newton / meter * angstrom * angstrom).rescale(eV)
    gradunits = (newton / meter * angstrom).rescale(eV/angstrom) / scale * 0.5
    for endpoint, vector in node.sc_bond_iter():
      vector = dot(cell, vector) + endpoint.pos - node.pos
      params = self[node.type, endpoint.type]
      bond_length = params[0]

      e0 = sum(vector*vector) * scale2 / bond_length - bond_length
      mult = params[1:] * bondparams * e0
      energy += e0 * (mult[0] + e0 * ( mult[1]                                    \
                 + e0 * ( mult[2] + (mult[3] + e0 * mult[4]) ) ))

      if stress is not None: 
        mult = params[1:] * gbondparams * e0
        e0grad = 2e0 * scale2 / bond_length                                    \
                 * (mult[0] + e0 * ( mult[1]                                   \
                   + e0 * ( mult[2] + (mult[3] + e0 * mult[4]) ) ))
        hold = e0grad * vector * gradunits
        node.center.gradient -= hold
        endpoint.center.gradient += hold
  
        matrix = outer(vector, vector)
        stress += e0grad * 0.5 * matrix * stressunits
                
    return energy * 3e0 / 8e0 * (newton/meter*angstrom*angstrom).rescale(eV)

  def _evaluate_angles(self, node, scale, cell, stress=None):
    """ Evaluates bond-angle and bond-bending modes. 

        .. note:: stress is not normalized by -1/volume yet.
    """
    from numpy import dot, sum, array, sqrt, outer
    from quantities import newton, meter, angstrom, eV
    bondparams = array([1e0, 5e-1/sqrt(3e0), 1e0/16e0, sqrt(3e0)/8e0, 1e0/640])
    gbondparams = array([ 0.75e0, 3e0*sqrt(3e0)/16.0, 3e0/32e0,
                          3e0*sqrt(3e0)/256e0, 0.003515625 ])
    scale2 = scale * scale

    energy = 0
    stressunits = (newton / meter * angstrom * angstrom).rescale(eV)
    gradunits = (newton / meter * angstrom).rescale(eV/angstrom) / scale * 0.5
    for (A,dA), (B, dB) in node.angle_iter():
      vA = dot(cell, dA) + A.pos - node.pos
      vB = dot(cell, dB) + B.pos - node.pos
      paramsA = self[node.type, A.type]
      paramsB = self[node.type, B.type]
      paramsAB = self[A.type, node.type, B.type]
      lengthA, lengthB = paramsA[0], paramsB[0]
      gamma, sigma =  paramsAB[0], paramsAB[1]
      mean_length = sqrt(lengthA * lengthB)

      e0 = sum(vA*vA) * scale2 / lengthA - lengthA                             \
            + sum(vB*vB) * scale2 / lengthB - lengthB
      e1 = dot(vA, vB) * scale2 / mean_length - mean_length * gamma 

      # bond-bending
      mult = paramsAB[2:] * bondparams * e1
      energy += e1 * (mult[0] + e1 * ( mult[1]                                 \
                 + e1 * ( mult[2] + (mult[3] + e1 * mult[4]) ) ))

      if stress is not None: 
        mult = paramsAB[2:] * gbondparams * e1
        e1grad = 2e0 * scale2 / mean_length                                    \
                 * (mult[0] + e1 * ( mult[1]                                   \
                   + e1 * ( mult[2] + (mult[3] + e1 * mult[4]) ) ))
        hold0 = e1grad * vA * gradunits
        hold1 = e1grad * vB * gradunits
        node.center.gradient -= hold0 + hold1
        A.center.gradient += hold1
        B.center.gradient += hold0
  
        matrix = outer(vA, vB)
        matrix += matrix.T
        stress += e1grad * 0.5 * matrix * stressunits

      # bond angle 
      energy += e0 * e1 * sigma
     
      if stress is not None: 
        hold0 = 1.5 * e1 * sigma / lengthA * scale2 * vA                       \
                + 0.75 * e0 * sigma / mean_length * scale2 * vB 
        hold0 = hold0 * gradunits
        hold1 = 1.5 * e1 * sigma / lengthB * scale2 * vB                       \
                + 0.75 * e0 * sigma / mean_length * scale2 * vA  
        hold1 = hold1 * gradunits


        node.center.gradient -= hold0 + hold1
        A.center.gradient += hold0
        B.center.gradient += hold1

        matrix = outer(vA, vB)
        matrix = 2e0 * e1 * (outer(vA, vA)/lengthA + outer(vB, vB)/lengthB)    \
                 + e0 / mean_length * (matrix + matrix.T)
        stress += matrix * 0.375 * sigma * scale2 * stressunits

    return energy * 3e0 / 8e0 * (newton/meter*angstrom*angstrom).rescale(eV)

