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

def inequivalent_sites(lattice):
  """ Returns a list containing only one site index per inequivalent sub-lattice. """
  from pylada.crystal import which_site

  result = set( i for i in range(len(lattice.sites)) ) 
  for i, site in enumerate(lattice.sites):
    if i not in result: continue
    for op in lattice.space_group:
      j = which_site( op(site.pos), lattice )
      if j != i and j in result: result.remove(j)
  return result

from numpy.linalg import norm
from pylada.crystal import A2BX4
from pylada.ce import create_clusters
from math import exp

# input
shell = 35

# create spinel lattice.
lattice = A2BX4.b5()
# occupations.
for site in lattice.sites:
  if "X" in site.type: continue
  site.type = ["A", "B"]
# recomputes space group for safety.
lattice.find_space_group()
lattice.set_as_crystal_lattice()

# creates random structure.
structure = A2BX4.b5I().to_structure() 
normal = A2BX4.b5().to_structure()
# for atom in structure.atoms:
#   atom.type = choice(lattice.sites[atom.site].type)


# list of inequivalent sites, keeping only those without "X"
ineqs = []
for i in  inequivalent_sites(lattice):
  if len(lattice.sites[i].type) > 1: ineqs.append(i)


# now creates multi-lattice clusters, along with index bookkeeping.
# first J0
clusters = create_clusters(lattice, nth_shell=0, order=0, site=0)
# then J1 for these sites.
for i in ineqs:
  clusters.extend(create_clusters(lattice, nth_shell=0, order=1, site=i))
# then pairs figures.
clusters = None
for i in ineqs:
  if clusters is None:
    clusters = create_clusters(lattice, nth_shell=shell, order=2, site=i)
  else: 
    clusters.extend(create_clusters(lattice, nth_shell=shell, order=2, site=i))
# sets interaction energies (eci)
# J0
clusters[0].eci = 0e0
# J1
for i in range(1, len(ineqs)+1): clusters[i].eci = 0e0
# J2
# clusters -> list of pair figures.
# clusters[i] -> a list of symmetrically *equivalent* figures.
# clusters[i][0] -> one single figure.
# clusters[i][0].origin -> origin of the figure: eg first spin.
# clusters[i][0][j] -> vector from origin to other spins:
#    eg clusters[i][0].origin.pos + clusters[i][0][0].pos = position of second spin
for i in range(len(ineqs)+1, len(clusters)): 
# print norm(clusters[i][0][0].pos)
  clusters[i].eci = exp(-5e0*norm(clusters[i][0][0].pos))

# now computes pis.
pis = clusters.pis(structure) - clusters.pis(normal)
# print pis
# now compute energies.
print shell, clusters(structure) - clusters(normal)
