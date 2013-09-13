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

""" Contains GA operators and evaluators for ESCAN properties. 

    At present operators are included only for elemental superlattices.
"""
__docformat__ = "restructuredtext en"
__all__ = ['elemental', 'CompareSuperCells']

class CompareSuperCells(object):
  """ A class to compare crystal supercells according to occupation. 
  
      Takes a supercell and a lattice on initialization:
      
      >>> comparison = CompareSuperCells(lattice, supercell)

      And then compares two structures as a functor call. It returns True if
      the two structures are comparable.

      >>> comparison(A, B)
  """
  def __init__(self, lattice=None, supercell=None, converter = None):
    """ Initializes comparison operator. 

        :Parameters:
          lattice : `pylada.crystal.Lattice`
            The backbone lattice.
          supercell : nump 3x3 array
            The supercell for which configurations are given.
          converter : None or callable 
            If a callable then will check input for attribute *integer_has_code*.
            If attribute does not exist, will convert input using this
            converter and obtain integer hash_code and resulting instance. If
            None, functor expects `pylada.crystal.Structure` instances which it
            will compare.
    """
    from ...crystal import smith_normal_transform

    self.lattice  = lattice
    """ Backbone lattice for which to compare supercells. """
    # check lattice is OK.
    self.check_lattice()
    
    self.supercell = supercell 
    """ Supercell for which to do comparisons. """
  
    self.converter = converter
    """ Converts individuals into structures. """

    try:
      from ...enumeration import count_flavors, create_flavorbase, create_transforms

      self.smith = smith_normal_transform(lattice.cell, supercell)
      """ Smith normal transform tuple """
      self.nflavors = count_flavors(lattice)
      """ Number of flavors. """
      self.nsites = len([0 for i in lattice.sites if len(i.type) > 1])
      """ Number of combinatorial sites. """
      self.card = int(self.smith[1][0]*self.smith[1][1]*self.smith[1][2]*self.nsites)
      """ Size of the configuration space. """
      self.flavorbase = create_flavorbase(self.card, self.nflavors)
      """ Basis of different flavors. """
      self.indices = []
      """ Indices to combinatorial sites. """
      i = 0
      for site in lattice.sites:
        if len(site.type) <= 1: self.indice.append(None)
        else: i += 1; self.indices.append(i-1)
      transforms = create_transforms(self.lattice)
      self.transforms = []
      """ Rotational transforms. """
      for transform in transforms:
        if not transform.invariant(self.supercell): continue
        transform.init(*self.smith)
        if not transform.is_trivial: self.transforms.append(transform)

      self.translations = Translation(self.smith[1], self.nsites)
      """ Translation transforms. """
      self._doenum = True
      """ True if should use enumeration package. """
    except: # supercell is likely too large or enumeration package not included.
      self._doenum = False


  
  def check_lattice(self):
    """ Raises an TypeError exception if the lattice is not usable. """
    sites = [site.type for site in self.lattice.sites if len(site.type) > 1]
    # should have a site with multiple occupation
    assert len(sites) > 0, TypeError("No combinator site found in lattice.")
    if len(sites) > 1:
      # not necessary if different wyckoff positions! Always the case?
      # still need same number of species though.
      for site in sites[1:]:
        assert len(set(site) ^ set(sites[0])) == 0,\
               TypeError( "Enumeration for combinatorial sites with "\
                          "different occupation not implemented. " )
        for a, b in zip(site, sites[0]):
          assert a == b, TypeError("Combinatorial sites have different order.")

  def structure_to_integer(self, structure):
    """ Given a structure, what is its integer hash-code? 
    
        The hash-code as defined here is not minimal, since the same
        structure but rotated/translated and otherwise transformed may result
        in a different integer number.
    """
    assert self._doenum, \
           RuntimeError( "Cannot create hash-code.\n"\
                         "Enumeration package absent, or structure too large.") 
    from numpy import abs, any
    from ...crystal import linear_smith_index

    # this is not an accurate definition. should use something like integer transform matrix...
    assert not any(abs(structure.cell - self.supercell) > 1e-12), \
           ValueError("Structure and supercell are not coherent.")

    # now creates integers from linear smith index
    result = 0
    for atom in structure.atoms: 

      assert atom.site >= 0 and atom.site <= len(self.lattice.sites),\
             ValueError("Atoms in structure are not correctly indexed by site type.")

      site = self.lattice.sites[atom.site]
      if len(site.type) < 2: continue

      lsi = linear_smith_index(self.smith, atom.pos - site.pos, self.indices[atom.site])
      assert atom.type in site.type,\
             ValueError( "Atomic type in structure not found in backbone lattice.\n{0} not in {1}"\
                         .format(atom.type, str(site.type)))
      result += self.flavorbase[lsi] * site.type.index(atom.type)

    return result

  def structure_hashcode(self, structure):
    """ Returns structure hashcode according to enumeration algorithms. """
    assert self._doenum, \
           RuntimeError( "Cannot create hash-code.\n"\
                         "Enumeration package absent, or structure too large.") 
         
    result = self.structure_to_integer(structure)
    x = int(result)

    # checks pure translation
    for translation in self.translations:
      t = translation(x, self.flavorbase)
      if t < result: result = t

    # checks rotation + pure translation
    for transform in self.transforms:
      t = transform(x, self.flavorbase)
      if t == x: continue
      if t < result: result = t
      for translation in self.translations:
        tt = translation(t, self.flavorbase)
        if tt < result: result = tt

    return result

  def __call__(self, a, b):
    """ Returns True if two structures are equivalent. """
    if not self._doenum: return not any(a.genes != b.genes)
    if hasattr(self.converter, "__call__"):
      if not hasattr(a, "integer_hashcode"): 
        a.integer_hashcode = self.structure_hashcode(self.converter(a.genes))
      if not hasattr(b, "integer_hashcode"): 
        b.integer_hashcode = self.structure_hashcode(self.converter(b.genes))
      return a.integer_hashcode == b.integer_hashcode

    else:
      return self.structure_hashcode(a) == self.structure_hashcode(b)

  def __getstate__(self):
    """ Saves state. """
    return self.lattice, self.supercell, self.converter
  def __setstate__(self, args):
    """ Resets state. """
    self.__init__(*args)

  @staticmethod
  def remove_hashcode(indiv):
    """ Removes hash-code from object. """
    if hasattr(indiv, "integer_hashcode"): del indiv.integer_hashcode
