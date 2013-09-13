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

""" Function object to convert to and from bitstrings and nanowires. """
__docformat__ = "restructuredtext en"
class Converter(object):
  """ Converts to and from bitstrings and nanowires. """

  def __init__( self, lattice, growth=(0,0,1), core_radius=2,\
                core_type='Si', types=['Si', 'Ge'], thickness=0.35,
                separation=3, passivant='Hg' ):
    """ Initializes a multi-shell wire.

        :Parameters:
          lattice 
            Lattice making up the backbone of the nanowire.
          growth : (0,0,1) or (1, 1, 0)
            Growth direction of the wire.
          core_radius : float
            radius of the core. 
          core_type : str
            Type of the core.
          types : list of str
            Types of atoms in the nanowire. It is used to determine what
            integer correspond to in a bitstring. 
          thickness : float
            Thickness of a shell. 
          separation : int
            Distance between the edge of a nanowire and the edge of the cell.
          passsivant : str
            Type of the passivant/barrier.
    """
    from copy import deepcopy

    self.lattice = lattice
    """ Zinc-Blende lattice for which the nanowire is constructed. """
    self.growth = growth
    """ Growth direction of the nanowire. """
    self.core_radius = core_radius
    """ Radius of the nanowire's core. """
    self.core_type = core_type
    """ Type of the core. """
    assert self.core_type in [u for a in lattice.sites for u in a.type],\
           ValueError("Could not find atomic-specie of core in lattice.")
    self.types = deepcopy(types)
    """ List of types which GA will optimize.

        Serves to determine what 0 and 1 correspond to in bitstring.
    """
    self.passivant = passivant
    """ Type of the passivant. """
    self.thickness = thickness
    """ Thickness of a shell in units of `Lattice.scale`. """
    self.separation = separation
    """ Number of shells between nanowire and end of cell. """

  def supercell(self, l_plane):
     """ Constructs supercell. """
     from numpy import array, all, abs, round, sqrt
     l_xy = int(round(l_plane/sqrt(2.0)+0.5))
     l_plane = int(round(l_plane+0.5))
     if all( abs(array(self.growth) - [0,0,1]) < 1e-12 ):
       cell = array([[l_plane, 0, 0],[0, l_plane, 0], [0, 0, 1.0]])
       return self.lattice.to_structure(cell)
     if all( abs(array(self.growth) - [1,1,0]) < 1e-12 ):
       cell = array([[-l_xy, 0, 0.5],[l_xy, 0, 0.5], [0, l_plane, 0]]) 
       return self.lattice.to_structure(cell)
     raise NotImplementedError("Unknown growth direction {0}.".format(self.growth))

  def all_shells(self, bitstring):
    """ Returns an array with all shells, including core. """
    from numpy import concatenate, zeros
    core = zeros(shape=(self.core_radius,), dtype=getattr(bitstring, 'dtype', 'int64'))
    core[:] = self.types.index(self.core_type)
    return concatenate((core, bitstring))

  def superstructure(self, shells):
    """ Returns an empty superstructure of the correct size.
    
        The atoms in the superstructure are correctly centered.
    """
    from numpy import dot, floor
    from numpy.linalg import inv
    # constructs the super-structure which will contain the nanowire
    structure = self.supercell( (len(shells) + self.separation) * self.thickness * 2.)
    
    # Centers atoms in the supercell of the superstructure.
    invcell = inv(structure.cell) 
    for atom in structure.atoms:
      # computes fractional coordinates and recenters atom to middle of cell.
      fractional  = dot(invcell, atom.pos) + [0.5, 0.5, 0] 
      # Move fractional coordinates to [0, 1] interval, taking into account numerical noise.
      fractional -= floor(fractional + 1e-12)
      assert all(fractional + 1e-12 >= 0e0) and all(fractional + 1e-12 < 1e0)
      # Move from fractional back to cartesian coordinates.
      atom.pos = dot(structure.cell, fractional)

    return structure
    
  def to_wire(self, bitstring):
    """ Returns a structure constructed from a bitstring. """
    from numpy import dot
    from ....crystal.iterator import shell as shell_iterator
    # creates an array of shell types, including nanowire core.
    shells = self.all_shells(bitstring)
    # create the superstructure which will contain the nanowire.
    nanowire = self.superstructure(shells)

    # finds the center of the cell (and nanowire).
    center = dot(nanowire.cell, [0.5,0.5,0])

    # now iterates over atoms in each shell.
    atomic_shells = shell_iterator(nanowire, center, self.growth, self.thickness)
    for shell_index, atoms in enumerate(atomic_shells):
      # finds the type, taking care whether we are inside or outside the nanowire.
      type = self.types[int(shells[shell_index])] if shell_index < len(shells) else self.passivant
      # sets the atoms in that layer to the correct type.
      for atom in atoms: atom.type = type

    return nanowire

  def to_bitstring(self, nanowire):
    """ Constructs a bitstring from a nanowire. """
    from numpy import array, dot
    from ....crystal.iterator import shell as shell_iterator

    result = []
    # finds the center of the cell (and nanowire).
    center = dot(nanowire.cell, [0.5,0.5,0])
    # loop over shells in nanowire.
    for layer in shell_iterator(nanowire, center, self.growth, self.thickness):
      # gets first atom in layer. 
      atom = layer.next() 
      # break from loop if we reached passivant.
      if atom.type == self.passivant: break 
      # add to list of types.
      result.append( self.types.index(atom.type) )

    return array(result[self.core_radius:])

  def __call__(self, object):
    """ Converts to and from structures and bitstring. """
    if hasattr(object, 'cell'): return self.to_bitstring(object)
    elif hasattr(object, 'genes'): return self.to_wire(object.genes)
    return self.to_wire(object)
 
  def __repr__(self):
    """ Returns string python code representing this object. """
    return "from {0.__class__.__module__} import {0.__class__.__name__}\n"\
           "converter = {0.__class__.__name__}( lattice, growth={1},\n"\
           "                 core_radius={0.core_radius}, core_type='{0.core_type}',\n"\
           "                 types={2}, thickness={0.thickness}, passivant={3}, \n"\
           "                 separation={0.separation} )"\
           .format(self, repr(self.growth), repr(self.types), repr(self.passivant))
