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

#! /usr/bin/python
def enum( _n, _lattice ):
  import numpy as np
  from pylada import enumeration
  from math import pow

  supercells = enumeration.find_all_cells(_lattice, _n)
  smiths = enumeration.create_smith_groups(_lattice, supercells)
  nflavors = enumeration.count_flavors(_lattice)
  nsites = len(_lattice.sites)
  transforms = enumeration.create_transforms(_lattice)
  for smith in smiths:
    card = int(smith.smith[0]*smith.smith[1]*smith.smith[2]*nsites)
    label_exchange=enumeration.LabelExchange( card, nflavors )
    flavorbase = enumeration.create_flavorbase(card, nflavors)
    translations = enumeration.Translation(smith.smith, nsites)
    database = enumeration.Database(card, nflavors)
    maxterm = 0
    for x in xrange(1, int(pow(nflavors, card))-1):
      if not database[x]: continue
      maxterm = x
      for labelperm in label_exchange:
        t = labelperm(x, flavorbase)
        if t > x: database[t] = False
      for translation in translations:
        t = translation(x, flavorbase)
        if t > x: database[t] = False
        elif t == x: 
          database[t] = False
          continue
        for labelperm in label_exchange:
          u = labelperm(t, flavorbase)
          if u > x: database[u] = False

    # checks supercell dependent transforms.
    for nsupercell, supercell in enumerate(smith.supercells):
      mine = []
      # creates list of transformation which leave the supercell invariant.
      cell = np.dot(_lattice.cell,  supercell.hermite)
      specialized = []
      for transform in transforms:
        if not transform.invariant(cell): continue
        transform.init(supercell.transform, smith.smith)
        if not transform.is_trivial:
          specialized.append( transform )

      specialized_database = enumeration.Database(database)
      for x in xrange(1, maxterm+1):
        if not database[x]: continue
        maxterm = x
        
        for transform in specialized:
          t = transform(x, flavorbase)
          if t == x: continue
          specialized_database[t] = False

          for labelperm in label_exchange:
            u = labelperm(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            
          for translation in translations:
            u = translation(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            for labelperm in label_exchange:
              v = labelperm(u, flavorbase)
              if v == x: continue
              specialized_database[v] = False
        if specialized_database[x]: yield x, smith, supercell, flavorbase

def list_all_structures( _n0, _n1 ):

  from pylada import crystal

  lattice = crystal.Structure().lattice

  result = []
  for n in range(_n0, _n1):
    oldsupercell = None
    for x, smith, supercell, flavorbase in enum(n, lattice):
      if oldsupercell is None or oldsupercell != supercell:
        result.append( (supercell, flavorbase, smith, []) )
        oldsupercell = supercell
      result[-1][3].append(x)
  return lattice, result

def choose_structures( _howmany, _min = 2, _max = 9 ):
  import random
  import numpy as np
  from pylada import crystal, enumeration

  lattice, str_list = list_all_structures( _min, _max )

  for i in  xrange(_howmany):
    nsupercell = random.randint(0, len(str_list)-1)
    while len(str_list[nsupercell][3]) == 0: 
      nsupercell = random.randint(0, len(str_list)-1)
    nstr = 0
    if len(str_list[nsupercell][3]) > 1: 
      nstr = random.randint(0, len(str_list[nsupercell][3])-1)

    supercell, flavorbase, smith = str_list[nsupercell][:3]
    x = str_list[nsupercell][3].pop(nstr)

    # creates structure
    structure = crystal.sStructure()
    structure.cell = np.dot(lattice.cell, supercell.hermite)
    crystal.fill_structure(structure)
    structure.scale = lattice.scale
    enumeration.as_structure(structure, x, flavorbase)
    yield crystal.Structure(structure)
    
def convert(_classes):
  from pylada import ce

  result = ce.ClusterClasses()
  for class_ in _classes:
    dummy = ce.Clusters()
    for mlcluster in class_:
      cluster = ce.Cluster()
      cluster.eci = class_.eci
      cluster.vectors.append(mlcluster.origin.pos)
      for i, vec in enumerate(mlcluster):
        cluster.vectors.append(vec.pos)
      dummy.append(cluster)
    result.append(dummy)
  return result

def main():
  from pylada import ce
  import boost.mpi as mpi
  from numpy import array
  from pylada.crystal.bravais import fcc

  functional = ce.Cubic()
  functional.set_mpi(mpi.world)
  functional.load( "input.xml" )
  lattice = fcc(); lattice.sites[0].type = "Au", "Pd"
  structure = lattice.to_structure( array([[0,1.0,1.5],[1,0,1.5],[1.0,1.0,0]]))
# print structure.lattice.space_group

  mlclasses = ce.MLClusterClasses("input.xml", False)
# mlclasses[0]
# print mlclasses
  classes = convert(mlclasses)
# print classes
# print
# print functional.clusters

# mlclusters = convert_clusters_to_mlclusters
  tests = [ ("000011100000", -74.377010),
            ("111000000101", -83.064845), 
            ("101000000011", -80.889828), 
            ("000111111000", -84.277151), 
            ("000111111001", -71.156430), 
            ("110000000111", -83.064845), 
            ("110000000101", -80.889828), 
            ("111000000111", -84.277151), 
            ("000111101000", -83.064845), 
            ("000110111000", -83.064845), 
            ("101000000111", -83.064845), 
            ("000110011000", -80.889828), 
            ("000011101000", -80.889828), 
            ("000111011000", -83.064845), 
            ("000011111000", -83.064845), 
            ("110000000011", -80.889828), 
            ("111000100111", -71.156430), 
            ("001000000000", -27.375282), 
            ("011000000111", -83.064845), 
            ("000010000000", -27.375282), 
            ("000000100000", -27.375282), 
            ("000000000001", -27.375282), 
            ("000100001000", -55.269849), 
            ("111010000111", -71.156430), 
            ("000000000100", -27.375282), 
            ("010000000000", -27.375282), 
            ("000111111100", -71.156430), 
            ("000001100000", -55.269849), 
            ("000011110000", -80.889828), 
            ("111001100111", -58.030174), 
            ("100000000010", -55.269849), 
            ("010000000100", -55.269849), 
            ("000001110000", -74.377010), 
            ("010000000001", -55.269849) ] 

  for test in tests:
    for i, atom in enumerate(structure.atoms):
      if test[0][i] == "0": atom.type = "Au"
      else:                 atom.type = "Pd"
    e = mlclasses(structure) 
    c = functional.chemical(structure) 
    if not e == 0e0:
      print "check %f, test %f, diff %f " \
            % ( c, e, c-e )
# for structure in choose_structures(15, 1, 12):
#   e = mlclasses(structure) 
#   c = functional.chemical(structure) 
#   if not e == 0e0:
#     print "check %f, test %f, diff %f " \
#           % ( e, c, fabs(c-e) )

#         % ( c, test[1], fabs(c-test[1]) )

if __name__ == "__main__":
  main()
