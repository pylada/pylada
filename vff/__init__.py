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

""" Valence Force Field for zinc-blende. 

 
    Implementation of the valence force field method for zinc-blende.
"""
__all__ = ['Node', 'Vff', 'read_input', 'exec_input', 'build_tree']
from .cppwrappers import Node
from .functional import Functional as Vff

def read_input(filepath="input.py", namespace=None):
  """ Specialized read_input function for vasp. 
  
      :Parameters: 
        filepath : str
          A path to the input file.
        namespace : dict
          Additional names to include in the local namespace when evaluating
          the input file.

      It add a few names to the input-file's namespace. 
  """
  from ..misc import read_input
  from . import specie
  from relax import Epitaxial, Relax

  # names we need to create input.
  input_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input': global_dict[k] = globals()[k]
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)

def exec_input( script, global_dict=None, local_dict=None,
                paths=None, name=None ):
  """ Specialized exec_input function for vasp. """
  from ..misc import exec_input

  # names we need to create input.
  if global_dict is None: global_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input': global_dict[k] = globals()[k]
  return exec_input(script, global_dict, local_dict, paths, name)

def build_tree(structure, overlap=1.2, **kwargs):
  """ Build first neighbor tree. """
  from numpy import array, argsort, dot, sum
  from numpy.linalg import inv
  from quantities import angstrom
  from ..crystal import periodic_dnc
  from ..error import ValueError

  if hasattr(overlap, 'rescale'):
    overlap = float(overlap.rescale(angstrom)) / structure.scale

  # net: map of atom id to the nodes in the the first neighbor net which
  #      wraps them.
  net = {}
  for i, atom in enumerate(structure): net[id(atom)] = Node(atom, i)
  
  # invcell: transform from cartesian to fractional coordinates.
  invcell = inv(structure.cell)
  # boxes: list of divide and conquer boxes.
  boxes = periodic_dnc(structure, overlap=overlap, **kwargs)
  # loop over boxes.
  for box in boxes:
    if len(box) < 5: raise ValueError('Box does not contain enough atoms.')
    # position: All positions in box as a numpy array, for faster treatment.
    positions = array([a.pos + b for a, b, c in box])
    # Loop over each node in the net in that particular box.
    for position, (atom, vector, inbox) in zip(positions, box):
      # If atom is outside box, loop.
      if not inbox: continue 
      # node: Node which wraps current atom
      node = net[id(atom)]
      if len(node) == 4: continue # already found all bonds.

      # distances: distance from this atom to all others in box.
      distances = sum((positions - position)**2, axis=1)
      # sortees: indices of four nearest neighbors in box. 
      sortees = argsort(distances)[1:5]
      
      # Now adds sortee as bonds.
      for i in sortees: 
        frac = dot(invcell, box[i][1] - vector)
        othernode = net[id(box[i][0])]
        if len(othernode) == 4: continue
        node.link(othernode, frac)
        if len(node) == 4: break
  # returns net as a list.
  return net.values()

from pylada import is_interactive
if is_interactive:
  def plot_tree(structure, overlap=0.5):
    """ Plots a tree using mayavi. """
    from enthought.mayavi import mlab
    from numpy import array, sum, any, dot, concatenate, argmin
    tree = build_tree(structure, overlap=overlap)
 
    # find original positions
    ids, originals = [], []
    for i, node in enumerate(tree):
      if id(node.center) not in ids: 
        ids.append(id(node.center))
        originals.append(node.pos)
    originals = array(originals)
    # now adds neighbors, including those outside the unit cell
    positions = originals.copy()
    others, connections = [], []
    for i, node in enumerate(tree): 
      for neighbor, vector in node:
        position = neighbor.pos + dot(structure.cell, vector)
        dist = sum((positions-position)**2, axis=1) 
        if any(dist < 1e-8): connections.append( (i, argmin(dist)) )
        else: 
          positions = concatenate((positions, position[None, :]))
          connections.append((i, positions.shape[0]-1))
           
        if any(sum((position[None,:]-originals)**2, axis=1) < 1e-8): continue
        if len(others)                                                           \
           and any(sum((position[None,:]-others)**2, axis=1) < 1e-8): continue
        others.append(position)
    others = array(others)
 
    a0, a1, a2 = structure.cell.T
    box = [ [0, 0, 0], a0, a0+a1, a1, [0,0,0], a2, a2+a0, a2+a1+a0, a2+a1, a2,  \
            a2 + a0, a0, a0+a1, a0+a1+a2, a1+a2, a1 ] 
    mlab.plot3d(*array(box).T, tube_radius=0.01, color=(0.8, 0.8, 0.8))
 
    mlab.points3d(*originals.T, color=(0.8, 0.5, 0.5), scale_factor=0.2)
    mlab.points3d(*others.T, color=(0.5, 0.8, 0.5), scale_factor=0.1)
    # add connections
    pts = mlab.points3d(*positions.T, scale_factor=0)
    pts.mlab_source.dataset.lines = array(connections)
    tube = mlab.pipeline.tube(pts, tube_radius=0.01)
    tube.filter.radius_factor = 1.
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    mlab.pipeline.surface(tube, color=(0.3, 0.3, 0.3))

  def plot_box(cell, tube_radius=0.01, color=None):
    """ Adds rectangular box to mayavi. """
    from enthought.mayavi import mlab 
    from numpy import array
    a0, a1, a2 = cell.T
    box = array([[0, 0, 0], a0, a0+a1, a1, a2, a2+a0, a2+a1+a0, a2+a1])
    connections = array([ [0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [4, 5], 
                          [5, 6], [6, 7], [7, 4], [1, 5], [7, 3], [6, 2] ])

    pts = mlab.points3d(*box.T, scale_factor=0)
    pts.mlab_source.dataset.lines = connections
  
    # Use a tube fiter to plot tubes on the link, varying the radius with the
    # scalar value
    tube = mlab.pipeline.tube(pts, tube_radius=tube_radius)
    tube.filter.radius_factor = 1.
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    if color is None: color = (0.8, 0.8, 0)
    mlab.pipeline.surface(tube, color=color)
del is_interactive



