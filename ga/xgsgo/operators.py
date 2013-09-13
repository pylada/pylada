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

def fractional_pos(structure, roll=True):
  """ Returns fractional coordinates of atoms in structure. """
  from random import random
  from numpy import dot, array
  from numpy.linalg import inv
  # random translation vector.
  trans = (random(), random(), random()) if roll else (0,0,0)
  result = [dot(inv(structure.cell), atom.pos) + trans for atom in structure]

  for pos in result: pos = mod(pos, [1,1,1])

  return array(result)

def cut_and_splice(s1, s2, roll=True):
  """ Cut-n-splice GSGO crossover operation

      Mating operation on two parent structures s1 and s2.
      Done on scaled atomic positions (cubic systems) by
      cutting them in half and mixing their upper and 
      lower parts.
  """
  from random import choice, random
  from numpy import dot, abs
  from numpy.linalg import det
  from pylada.crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos1 = zip([atom.type for atom in s1],fractional_pos(s1, roll=True))
  sc_pos2 = zip([atom.type for atom in s2],fractional_pos(s2, roll=True))

  result = Structure(s1.cell, scale=s1.scale)

  # choose random positions of split-plane
  xsep = 0.5 - (random() * 0.45 + 0.15)
  # choose direction of split-plane randomly from cell-vectors.
  direction = choice(range(3))

  for type, pos in sc_pos1:
    if pos[direction] >= xsep: result.add_atom(*dot(result.cell, pos), type=type)

  for type, pos in sc_pos2:
    if pos[direction] < xsep: result.add_atom(*dot(result.cell, pos), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)

  return result


def mix_poscars(s1,s2, roll=True):
  """ Crossover operations where atoms are from one and cell from other.

      Mating operation on two parent structures s1 and s2.
      Done on scaled atomic positions (cubic systems) by
      interchanging their cells and scaled positions.
      Returns two offspring structures each with the same 
      number of atoms as the parent from which the atoms 
      are inhereted.
  """
  from random import choice
  from numpy import dot
  from numpy.linalg import det
  from ..crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos2 = zip([atom.type for atom in s2],fractional_pos(s2, roll))

  # cell from s1
  result = Structure(s1.cell, scal=s1.scale)

  # atoms from s2
  for type, pos in sc_pos2:
    result.add_atom(*dot(result.cell, pos), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)
  return result

def mix_atoms(s1, s2, roll=True):
  """ Randomly mix cell and atoms from parents.

      Mating operation on two parent structures s1 and s2.
      Done by mixing randomly atoms from s1 and s2 into the 
      cell coming from one of them. Returns two offspring 
      structures.
  """
  from random import choice
  from itertools import chain
  from numpy.linalg import det
  from numpy import dot, abs
  from ..crystal import Structure

  # swap structures randomly.
  if choice([True, False]): s1, s2 = s2, s1

  # chem. symbols and scaled positions of the two parents
  sc_pos1 = zip([atom.type for atom in s1], fractional_pos(s1, roll))
  sc_pos2 = zip([atom.type for atom in s2], fractional_pos(s2, roll))

  result = Structure(s1.cell)

  for pos, type in chain(sc_pos1, sc_pos2):
    if choice([True, False]):
      result.add_atom(*dot(result.cell, pot), type=type)

  result.scale =  s1.scale * (float(len(result)) / float(len(s1)))**(1./3.)

  return result

def jiggle_structure(structure):
  """ Mutates by jiggling structure. """
  from copy import deepcopy
  from numpy import dot, min
  from numpy.linalg import norm, det
  from numpy.random import rand
  result = deepcopy(structure)
  omega = rand(3,3) * 2e0 - 1e0
  omega += omega.T
  newcell = result.cell + dot(omega, result.cell)

  fractionals = fractional_pos(result)
  mindists = [min([norm(a.pos - atom.pos) for a in result]) for atom in result]
  for atom, frac, mindist  in zip(result, fractionals, mindists):
    atom.pos = dot(newcell, frac)  + mindist * 0.25 * (rand(3) * 2e0 - 1e0)
  result.cell = newcell

  result.scale = structure.scale * (abs(det(structure.cell)) / abs(det(result.cell))\
                  * float(len(result)) / float(len(structure)))**(1./3.)
  return result
