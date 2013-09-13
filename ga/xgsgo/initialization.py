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

def random_cell(angle_range=None, cubic=False, **kwargs):
  """ Generates random cell. 
      
      In the algorithm below, a,b,c are the lengths of the lattice vectors
      alpha, beta, gamma are the angles betweeen them - the standard convetnion
      is respected (see `here <http://en.wikipedia.org/wiki/Crystal_structure>`_)
      lengths are between 0.6 and 1.4 in units of scale angles in the
      angle_range range

      :param float angle_range: 
           Minimum and maximum angles between lattice vectors.
           Defaults to 60.0 and 140.0 degrees.
      :param bool cubic: 
           Whether the cell should be cubic or not.
      :param kwargs:
           Passed on to :py:func:`~pylada.math.gruber`

      :returns: A random crystal cell.
  """
  from random import random
  from numpy import pi, cos, sin, array, zeros, abs, sqrt
  from numpy.linalg import det
  from ...math import gruber

  if angle_range == None: angle_range = [60.,140.]
  cell = zeros([3,3])
  while abs(det(cell)) < 1e-6:

    if cubic:
      a = b = c = 1.
      alpha = beta = gamma = pi/2.
    else:
      a = 0.8*random()+0.6 
      b = 0.8*random()+0.6
      c = 0.8*random()+0.6
        
      alpha = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.
      beta  = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.
      gamma = ((angle_range[1]-angle_range[0])*random() + angle_range[0])*pi/180.

    a1 = a*array([1.,0.,0.])
    a2 = b*array([cos(gamma),sin(gamma),0.])
        
    c1 = c*cos(beta)
    c2 = c/sin(gamma)*(-cos(beta)*cos(gamma) + cos(alpha))
        
    a3 = array([c1, c2, sqrt(abs(c**2-(c1**2+c2**2)))])
        
    cell = array([a1,a2,a3]).T
    if det(cell) < 0e0: cell = array([a2, a1, a3]).T
  return gruber(cell.copy(), **kwargs)

def random_structure( nbatoms, min_distance=1.3, atomic_scale=None,
                      specie='anion', **kwargs): 
  """ Generates random crystal structure.
      
      Creates a random cell using :py:meth:`random_cell`, then places atoms
      randomly inside the cell.
      contraints on the range of angles between lattice 
      vectors and the minimal distance distances between 
      two atoms.
      
      :param int nbatoms: 
         Number of atoms in the structure.
      :param real min_distance: 
         Minimal distance between two atoms. If no units, defaults to angstrom.
      :param real atomic_scale:
         Used in probability grid when an atom is placed. A guassian with this
         half-width is created which makes placing an atom here less probable.
         If None, then we use min_distance.
      :param **kwargs: 
         Extra parameters are passed on to :py:meth:`random_cell`.

      :returns: A random crystal structure.
  """
  from random import choice
  from numpy import dot, array, arange, ones, argwhere, exp, sum, argmin
  from itertools import product, chain
  from quantities import angstrom
  from ...crystal import Structure, neighbors

  # creates structure with random cell.
  result = Structure(random_cell(**kwargs))
  # sets scale using min_distance as a first guess.
  if hasattr(min_distance, 'rescale'):
    min_distance = float(min_distance.rescale(angstrom))
  minvol = 4.0/3.0*min_distance**3 * nbatoms
  result.scale = (2*minvol/result.volume)**(1.0/3.0)
  # set atomic_scale to result.scale units.
  if hasattr(atomic_scale, 'rescale'): 
    atomic_scale = float(atomic_scale.rescale(angstrom))
  if atomic_scale is None: atomic_scale = min_distance
  atomic_scale /= result.scale

  # Creates a grid where atoms will be placed.
  N = max(7, (3.0*nbatoms**(1./3.)))
  grid = array(list(product( arange(0, 1., 1.0/float(N)), 
                             arange(0, 1., 1.0/float(N)), 
                             arange(0, 1., 1.0/float(N)) )))
  grid = dot(result.cell, grid.T).T
  
  # Setting the initial prob. distribution
  probability  = ones(len(grid), dtype="float64")
  
  # periodic images. used to construct periodict probability density where
  # atoms shouldn't be.
  periodic_images = [ dot(result.cell, [i,j,k]) for i in arange(-1,2)\
                      for j in arange(-1,2) for k in arange(-1,2) ]
  periodic_images = array(periodic_images)

  # places first atom at origin.
  result.add_atom(0, 0, 0, type=specie)
  # make neighbors less probable
  distances = grid[:,None, :] - periodic_images[None, :, :] 
  probability -= sum(exp(-sum(distances**2, axis=2)/atomic_scale**2), axis=1)

  # now adds each atom in turn. 
  for i in xrange(1, nbatoms):
    # create criteria for placing atoms.
    maxi, mini = max(probability), min(probability)
    criteria = 0.90 * (maxi - mini) + mini
    # Find points on the anion grid that have the desired probability
    empty_space=argwhere(probability >= criteria) 
    # choose position randomly
    pom = choice(empty_space.flatten())
    # add atom to structure
    result.add_atom(*grid[pom].flatten(), type=specie)
    # make neighbors less probable
    distances = (grid- grid[pom])[:,None, :] - periodic_images[None, :, :] 
    probability -= sum(exp(-sum(distances**2, axis=2)/atomic_scale**2), axis=1)

  # Set the scale so that the min_distance is fulfilled
  mindist = [ (vector, distance) for atom in result \
              for neighbor, vector, distance in neighbors(result, 1, atom) ]
  mindist = mindist[argmin([d[1] for d in mindist])]
  if mindist[1] * result.scale < min_distance:
    result.scale *= 1.01*min_distance / result.scale / mindist[1] 

  return result

def populate_anion_cation(structure, species, anions):
  """ Populates structure with anions and cations. 
  
      Tries to place cations in such a way as to break as many anion-anion
      bonds. Then places the anion and cation species randomly within each
      sublattice.
  
      :param structure:
        :py:class:`~pylada.crystal.Structure` to modify. All atomic types are
        replaced.
      :param dict species: 
        Dictionary where keys correspond to atomic species, and values to the
        stoichiometry.
      :param dict anions:
        List of anionic species. Cations are the keys in ``species`` not listed
        in ``anions``.
  """
  from random import shuffle, choice
  from numpy import array, sum, argmax
  from ..crystal import coordination_shells

  # first creates list of cations and anions.
  types = anions
  anions, cations = [], []
  for key, value in species.iteritems():
    if key not in types: 
      cations.extend([key]*value)
    else: anions.extend([key]*value)
  if len(anions) + len(cations) != len(structure): 
    raise ValueError( "Number of species and atoms do not match.\n"\
                      "There are {0} anions and {1} cations, but "\
                      "{2} atoms."\
                      .format(len(anions), len(cations), len(structure)) )
    
  # makes sure all atoms are called anions.
  if len(anions) == 0: # unless no anions are requested.
    for a in structure: a.type = 'cation'
  else: 
    for atom in structure: atom.type = 'anions'
  # if cations exist, place them so as to minimize the number of cation-cation
  # and anion-anion bonds.
  if len(cations) != 0:
    # loop over shells, to minimize same-same bonds.
    shells = [ coordination_shells(structure, 1, atom, 0.5/structure.scale)[0] \
               for atom in structure ]
    # Each iteration places another cation.
    for i in xrange(len(cations)):
      # environment holds the same-same bond score.
      # since that score changes each time a cation is placed, it is
      # recomputed.
      environment = []
      for shell, atom in zip(shells, structure):
        if atom.type == 'cation': x = -1 # avoids double placement.
        else: 
          # score if function of bond-length.
          x = [(1.0/d if atom.type == n.type else 0.) for n, v, d in shell]
        environment.append(sum(x))
      environment = array(environment)
      # finds all highest scoring environments, and chooses a highest-scorer
      # randomly.
      indices = environment.argsort(axis=0)
      u = argmax(environment[indices])
      u = choice(indices[u:])
      structure[u].type = 'cation'

  # now replaces 'anion'/'cation' with actual species.
  shuffle(anions)
  shuffle(cations)
  for atom in structure:
    if atom.type == 'cation':
          atom.type = cations.pop()
    else: atom.type = anions.pop()

def populate(structure, species):
  """ Randomly populates structure with species. """
  from random import shuffle
  types = []
  for key, value in species.iteritems(): types.extend([key]*value)
  for atom, type in zip(structure, types): atom.type = type

def taboo( structure, max_atoms=-1, min_distance=1.3, \
           angle_range=None, same_first_neigh=2, verbose=True ):
  """ Checkes geometry of the structure and returnes 

      True if the structure satisfyies a set of conditions
      otherwise returns False
      
      :param structure:
          :py:class:`~pylada.crystal.Structure` to check.
      :param int max_atoms: 
          Maximal allowed number of atoms
      :param float min_distance: 
          Minimal first shell distance in angstroms
      :param angle_range: 
          Range of allowed angles between lattice vectors
      :type angle_range: 2x1 float list
      :param int same_first_neigh: 
          Maximal number of atoms that have more than two 
          neighbors of the same type whithin 1.5*min_distance.
          If strictly negative, then does not check.
      :param bool verbose: 
          Whether to print messages
  """
  from numpy import pi, arccos, dot, sqrt
  from ...crystal import neighbors
  
  # check number of atoms.
  if max_atoms > 0 and len(structure) > max_atoms: 
    if verbose: print "more than ",max_atoms," atoms"
    return True

  # chekc first neighbor distances.
  distances = []
  for atom in structure:
    for n, v, d in neighbors(structure,1,atom.pos):
      distances.append(d)

  if structure.scale*min(distances) < min_distance: 
    if verbose: print "distances shorter than",min_distance," AA"
    return True

  # check that angles are within specified range.
  if angle_range is None: angle_range = (45.,160.)
  cell = structure.cell
  alpha = 180./pi*arccos( dot(cell[:,0],cell[:,2]) /\
             sqrt(abs(dot(cell[:,0],cell[:,0])*dot(cell[:,2],cell[:,2])) ))
  beta  = 180./pi*arccos( dot(cell[:,1],cell[:,2]) /\
             sqrt(abs(dot(cell[:,1],cell[:,1])*dot(cell[:,2],cell[:,2])) ))
  gamma = 180./pi*arccos( dot(cell[:,0],cell[:,1])/\
             sqrt(abs(dot(cell[:,0],cell[:,0])*dot(cell[:,1],cell[:,1])) ))
  
  if not angle_range[0] <= alpha <= angle_range[1]:
    if verbose: print "alpha =",alpha," degrees"
    return True

  if not angle_range[0] <= beta  <= angle_range[1]:
    if verbose: print "beta =",beta," degrees"
    return True

  if not angle_range[0] <= gamma <= angle_range[1]:
    if verbose: print "gamma =",gamma," degrees"
    return True


  # check neighbor types.
  if same_first_neigh >= 0:
    check = []
    for atom in structure:
      pom = 0
      for n, v, d in neighbors(structure,20,atom.pos):
        if structure.scale*d <= 1.5*min_distance and n.type==atom.type:
          pom = pom + 1
      
      check.append(pom)

    how_many = len([x for x in check if x>2])

    if how_many > same_first_neigh: 
      if verbose: print how_many," atoms with more than 2 close neighbors of the same type"
      return True

  # all done
  return False


