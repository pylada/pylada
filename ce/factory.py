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

def cluster_factory(lattice, J0=False, J1=False, **mb):
  """ Returns class of cluster classes. 

      Creates a list of clusters using a simple input. 

      :param lattice:
          The backbone lattice on which the Ising model is created.
      :param bool J0:
          If True, adds zero order term.
      :param bool J1: 
          If True, adds on-site terms.
      :param mb: 
          All other many body terms, where the keys should be "B2", "B3", "Bn",
          where n is the order of the many-body interaction, and the values is
          the number of shells to look for.
  """
  from operator import itemgetter
  from numpy import dot, concatenate, sum
  from ._single_site_factory import check_keywords, factory as ss_factory
  from .cluster import Cluster, spin
  from ..crystal import which_site, space_group, coordination_shells

  # checks that keywords are well formed.
  check_keywords(J0, J1, **mb)

  if len(lattice) == 1: return ss_factory(lattice, J0, J1, **mb)

  # add space group if it does not already exist.
  if not hasattr(lattice, 'spacegroup'): 
    lattice.spacegroup = space_group(lattice)

  # computes equivalent sites and add stuff to lattice for later
  nonequiv = set( i for i in range(len(lattice)) ) 
  for i, site in enumerate(lattice):
    if i not in nonequiv: continue
    for op in lattice.spacegroup:
      j = which_site( dot(op[:3], site.pos) + op[3], lattice )
      if j != i and j in nonequiv: nonequiv.remove(j)
  lattice = lattice.copy()
  for i, site in enumerate(lattice):
    site.asymmetric = i in nonequiv
    site.index = i

  result = []
  # First creates J0
  if J0:
    result.append(Cluster(lattice))
  # then creates J1 for these sites.
  if J1: 
    for i in nonequiv: 
      site = lattice[i]
      if not hasattr(site.type, '__iter__'): continue
      result.append(Cluster(lattice))
      result[-1].spins = spin([0,0,0], site.index)

  if len(mb) == 0: return result
  # figure out max number of shells.
  maxshells = max(mb.itervalues())
  # creates many bodies.
  for i in nonequiv:
    site = lattice[i]
    if not hasattr(site.type, '__iter__'): continue
    # creates list of neighbors with indices into shells.
    allshells = coordination_shells(lattice, maxshells, site.pos)
    shells = []
    shellindices = []
    for shell in allshells:
      for atom, vector, distance in shell:
        if not hasattr(atom.type, '__iter__'): continue
        if not hasattr(atom.type, '__len__'): continue
        if len(atom.type) < 2: continue
        pos = vector.copy()
        if atom is not site: pos += site.pos - atom.pos
        shells.append(spin(pos, atom.index))
      shellindices.append(len(shells))
    shells = concatenate(shells)
    norms = sum(shells['position']* shells['position'], axis=1)
    sortee = [ (s, l, p[0], p[1], p[2])                                        \
               for s, l, p in zip( norms, shells['sublattice'], 
                                   shells['position']) ]
    sortee = [u[0] for u in sorted(enumerate(sortee), key=itemgetter(1))]
    shells = shells[sortee].copy()
    # loop over manybody figures and creates them.
    for key in sorted(mb.keys()):
      if mb[key] <= 0: continue
      order = int(key[1:])
      dummy = _create_clusters( lattice, 
                                neighbors=shells[:shellindices[mb[key]-1]],
                                site=site,
                                order=order )
      result.extend(dummy)
  return result

def _create_clusters(lattice, neighbors, site, order):
  """ Creates a set of clusters of given order up to given shell. """
  from numpy import any, concatenate, not_equal, array
  from .cluster import Cluster, spin
  from .cppwrappers import ProductILJ
  knowns = set()
  results = []
  # loop over figures made from those neighbors.
  # each figure is represented by an ordered set of indices in the neighbor
  # list.  This makes it easy to recognize which figures are already present in
  # the list of clusters.
  for indices in ProductILJ(range(len(neighbors)), order-1):
    # make sure this figure was not previously found by symmetry.
    if indices in knowns: continue
    # now create cluster
    a = Cluster(lattice)
    if order == 2:
      a.spins = array( [spin([0,0,0], site.index), neighbors[indices]],
                       dtype=neighbors.dtype )
    else: 
      a.spins = concatenate(( spin([0,0,0], site.index),
                              neighbors[list(indices)]), axis=0 )
    results.append(a)
    
    # now we add to the list of known figures those which are symmetrically
    # equivalent to this one.
    a._create_symmetrized()
    for u in a._symmetrized: 
      if u[0]['sublattice'] != site.index: continue
      indices = []
      for vector in u[1:]:
        for i, neigh in enumerate(neighbors):
          if vector['sublattice'] != neigh['sublattice']: continue
          if all(abs(vector['position'] - neigh['position']) < 1e-8):
            indices.append(i)
            break;
      assert len(indices) == order - 1, (u[1:], neighbors)
      indices = tuple(sorted(indices))
      if order == 2 or all(u < v for u, v in zip(indices[:-1], indices[1:])):
        knowns.add(indices)
  return results
