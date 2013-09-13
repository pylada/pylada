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

__docformat__ = "restructuredtext en"
__all__ = ['BasisSet', 'Shell']
from collections import MutableMapping
from quantities import UnitQuantity, angstrom
from .input import AttrBlock
from ..tools.input import VariableListKeyword, BaseKeyword

crystal_bohr = UnitQuantity('crystal_bohr', 0.5291772083*angstrom, symbol='bohr')
""" Bohr radius as defined by CRYSTAL. """
crystal_invbohr2 = UnitQuantity( 'crystal_invbohr2',
                                 1e0/crystal_bohr/crystal_bohr, symbol='bohr^-2' )

def specie_name(specie):
  """ Translate from CRYSTAL to sensible specie name. 

      CRYSTAL uses a fairly obscure numbering system to map species to
      basis-sets. This function translates that mapping into a more explicit
      naming scheme.
  """
  from ..periodic_table import find as find_specie
  try: n = int(specie)
  except: pass
  else: specie = n
  if isinstance(specie, str): 
    specie = specie.rstrip().lstrip()
    try: n = find_specie(name=specie)
    except: pass
    else: specie = n.symbol
  elif specie < 100:
    try: n = find_specie(atomic_number=specie)
    except: pass
    else: specie = n.symbol
  return specie

class Shell(object):
  """ Defines a general gaussian basis set for a specific orbital shell. 
  
      This class encodes a set of basis functions as follows.

      .. code-block:: python
        
         Shell('s', a0=[16120.0,  0.001959],
                    a1=[ 2426.0,  0.01493], 
                    a2=[  553.9,  0.07285],
                    a3=[  156.3,  0.2461], 
                    a4=[   50.07, 0.4859],
                    a5=[   17.02, 0.325], charge=0.0) 
      
      :param str type: 
         The first argument can be 's', 'p', 'd', or 'sp' and defines the
         channel for the set of orbitals. 
      
      :param float charge:
          Number of electrons assigned to this orbital at the start of the
          calculation. It should not be greater than the number of electrons
          that particular orbital channel can accomodate. For 's' it defaults
          to 2.0, for 'p' to 6, for 'd' to 10, and for 'sp' to 8.

      :param a0:
          First set of paramaters. For 's', 'p', 'd' channels, this should be a
          sequence of two floats, where the first float is the exponent of the
          gaussian and the second the contraction coefficient.
          In the case of 'sp' orbitals, it should be a sequence of three
          floats, where the first is the common exponent, and the second and
          third are the contraction coefficients for 's' and 'p' orbitals
          respectively.

      :param aN: 
          N-th set of orbital parameters.
  """
  def __init__(self, type='s', charge=None, **kwargs):
    """ Creates a gaussian basis set. 
    
        :params str type:
          Denotes the orbital type of this shell.
        :param float charge:
          Denotes the formal electronic charge on this shell. 
          Defaults(None) to a full shell.
    """
    from operator import itemgetter
    super(Shell, self).__init__()
    self.type = type
    """ Orbital type of this shell """
    self.charge = charge
    """ Nominal charge of the shell """
    self.functions = []
    """ Radial functions description """
    # collect basis function keywords.
    keys = []
    for key, value in kwargs.iteritems():
      if key[0] != 'a':
        raise ValueError( 'Unexpected keyword {0!r} in Shell.__init__(...)'    \
                          .format(key) )
      try: int(key[1:])
      except:  
        raise ValueError( 'Unexpected keyword {0!r} in Shell.__init__(...)'    \
                          .format(key) )
      try: value = list(value)
      except: 
        raise ValueError( 'Unexpected keyword value for {0!r} '                \
                          'in Shell.__init__(...): {1!r}'.format(key, value) )
      if len(value) != (3 if self.type == 'sp' else 2):
        raise ValueError( 'Unexpected sequence length for {0!r} '              \
                          'in Shell.__init__(...): {1!r}'.format(key, value) )
      keys.append((key, int(key[1:]), value))
    for key, i, value in sorted(keys, key=itemgetter(1)):
      self.append(*value)
  @property
  def type(self):
    """ Orbital shell of this function. """
    return self._type
  @type.setter
  def type(self, value):
    from ..error import ValueError
    if isinstance(value, int):
      if value == 0:   value = 's'
      elif value == 1: value = 'sp'
      elif value == 2: value = 'p'
      elif value == 3: value = 'd'
      elif value == 4: value = 'f'
      else: value = int(value)
    elif value not in ['s', 'p', 'd', 'f', 'sp']:
      try: value = int(value)
      except: raise ValueError('Wrong input to type ({0}).'.format(value))
      else:
        self.type = value
        return 
    self._type = value
  @property
  def type_as_int(self):
    """ Orbital type as CRYSTAL input. """
    if self._type == 's': return 0
    elif self._type == 'sp': return 1
    elif self._type == 'p': return 2
    elif self._type == 'd': return 3
    elif self._type == 'f': return 4

  @property
  def charge(self):
    """ Formal electronic charge attributed to this shell. """
    return self._charge
  @charge.setter
  def charge(self, value):
    if value is None or value < 0:
      if self._type == 's': value = 2
      elif self._type == 'sp': value = 8
      elif self._type == 'p': value = 6
      elif self._type == 'd': value = 10
      elif self._type == 'f': value = 18
    self._charge = float(value)

  def append(self, alpha, coef1, coef2=None): 
    """ Appends a gaussian to the basis set. """
    from ..error import input, TypeError
    if len(self.functions) >= 10: 
      raise input('Cannot input further coefficients.') 
    if self._type == 'd' and len(self.functions) >= 6:
      raise input('Cannot input further coefficients.') 
    if hasattr(alpha, 'rescale'): alpha = alpha.rescale(crystal_invbohr2)
    else: alpha = float(alpha) * crystal_invbohr2
    if self._type == 'sp': 
      if coef2 is None:
        raise TypeError('Expected second coeff for sp contraction.')
      self.functions.append([alpha, float(coef1), float(coef2)])
    else:
      if coef2 is not None:
        raise TypeError('Cannot use second coeff for in non-sp contraction.')
      self.functions.append([alpha, float(coef1)])
    return self

  @property
  def raw(self):
    """ Prints CRYSTAL input. """
    result = '0 {0.type_as_int} {1} {0.charge} 1\n'                            \
             .format(self, len(self.functions))
    for function in self.functions:
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbohr2).magnitude)
      if len(function) == 2:
        result += '{0:<20.12f}    {1: 20.12f}\n'.format(alpha, coef1)
      else:
        result += '{0:<20.12f}    {1:> 20.12f}  {2[2]:> 20.12f}\n'               \
                  .format(alpha, coef1, function)
    return result

  @raw.setter
  def raw(self, value):
    """ Reads input lines. """
    lines = value.split('\n')
    dummy, self.type, n, self.charge, dummy = lines.pop(0).split()
    self.functions = []
    for i in xrange(int(n)):
      line = lines.pop(0).split()
      self.append(*line[:3 if self._type == 'sp' else 2]) 

  def __repr__(self, indent=''):
    """ Represention of this instance. """
    args = ['{0.type!r}'.format(self)]
    if self.charge != {'s': 2, 'p': 6, 'sp': 8, 'd': 10, 'f': 18}[self.type]:
      args.append('{0.charge!r}'.format(self))
    for i, function in enumerate(self.functions): 
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbohr2).magnitude)
      vals = [alpha, coef1]
      if len(function) == 3: vals.append(function[2])
      args.append('a{0}={1!r}'.format(i, vals))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __len__(self): return len(self.functions)

class Ghosts(VariableListKeyword):
  """ Removes atoms but keeps basis functions. 

      This keyword should be in the basis subsection of the input.
      In practice, it can be set as follows::

        functional.basis.ghosts = [1, 3, 5], False

      The second item means that symmetries are *kept*.
      This results in the following crystal input:

        | KEEPSYMM
        | GHOSTS
        | 3
        | 1 3 5

      The first item in the tuple is a sequence of atomic lables (e.g.
      integers). The second is whether or not to break symmetries. It is True
      by default. In the above, symmetries are kept. Hence atoms equivalent
      those explicitely given are also removed.

      The two attributes can also be obtained and set using
      :py:attr:`~Ghosts.value` and :py:attr:`~Ghosts.breaksym':

        >>> functional.basis.ghosts.breaksym
        True
        >>> functional.basis.ghosts.value
        [1, 3, 5], 

  """
  type = int
  """ Type of the labels. """
  keyword = 'ghosts'
  """ CRYSTAL keyword """
  def __init__(self, value=None, **kwargs):
    """ Creates Ghost keyword. 
    
        :param [int] value:
          List of atom lables referencing atoms from which to make ghosts.
        :param bool keepsym: 
          Whether to keep symmetries or not. Defaults to True. If symmetries
          are kept, then all symmetrically equivalent atoms are ghosted. 
        :param bool breaksym:
          Whether to break symmetries or not. Defaults to False.
          Only one of breaksymm needs be specified. If both are, then they
          should be consistent.
    """
    super(Ghosts, self).__init__(value=value)
    self.breaksym = kwargs.get('breaksym', not kwargs.get('keepsym', True))
    """ Whether or not to break symmetries.
    
        Defaults to False. If symmetries are not broken, then all equivalent
        atoms are removed.
    """
  @property
  def keepsym(self):
    """ Not an alias for breaksym. """
    return not self.breaksym
  @keepsym.setter
  def keepsym(self, value):
    self.breaksym = not value
  def __repr__(self): 
    """ Dumps representation to string. """
    args = []
    if self.value is not None: args.append(repr(self.value))
    if self.breaksym == True: args.append("breaksym=True")
    return "{0.__class__.__name__}(".format(self) + ', '.join(args) + ')'
  def output_map(self, **kwargs):
    """ Prints CRYSTAL input """
    from ..tools.input import Tree
    if self.value is None: return None
    if len(self.value) == 0: return None
    if kwargs.get('breaksym', True) != self.breaksym:
      results = Tree() 
      results['keepsymm' if self.keepsym else 'breaksym'] = True
      results[self.keyword] = self.raw
      return results
    return {self.keyword: self.raw}
  def read_input(self, value, owner=None, **kwargs): 
    """ Reads from input. """
    self.raw = value
    self.breaksym = kwargs.get('breaksym', False)
  def __get__(self, instance, owner=None):
    return self
  def __set__(self, instance, value): 
    """ Setting Ghosts made easy. """
    if hasattr(value, '__iter__'): 
      self.value = list(value)
      return
    elif value is True or value is False: 
      self.breaksym = value
    elif hasattr(value, '__len__'):
      if len(value) == 2 and value[1] is True or value[1] is False:
        self.value, self.breaksym = value
      else: self.value = value
    else: raise ValueError( 'Ghosts expects a value of the form '              \
                            '``[integers], True or False``.')


class Chemod(BaseKeyword, MutableMapping):
  keyword = 'chemod'
  def __init__(self, **kwargs):
    MutableMapping.__init__(self)
    self.modifications = {}
    """ Holds chemical modifications. """
    self.breaksym = kwargs.get('breaksym', not kwargs.get('keepsymm', False))
    """ Whether to break symmetries or not. """
  @property
  def keepsymm(self):
    """ Alias to not breaksym. """
    return self.breaksym == False
  @keepsymm.setter
  def keepsym(self, value): self.breaksym = value == False

  def __getitem__(self, index): return self.modifications[index]
  def __setitem__(self, index, value): self.modifications[index] = value
  def __delitem__(self, index): self.modifications.__delitem__(index)
  def __len__(self): return len(self.modifications)
  def __iter__(self): return self.modifications.__iter__()
  def __contains__(self, name): return self.modifications.__contains__(name)
  def _symequivs(self, structure):
    """ Symetrically equivalent sites. """
    from numpy import dot
    from numpy.linalg import inv
    from ..crystal import which_site
    from ..error import internal
    result = {}
    sites  = structure.eval()
    symops = structure.symmetry_operators
    invcell = inv(sites.cell)

    # loop over atoms of interest
    for label in self.modifications:
      atom = sites[label-1]
      if atom.label != label: 
        raise ValueError( 'Incorrect atomic label for chemod ({0}).'           \
                          .format(label) )
      intermediate = result.get(label, set([label]))
      # look for symmetrically equivalents.
      for op in symops:
        newpos = dot(op[:3], atom.pos) + op[3]
        index = which_site(newpos, sites, invcell)
        if index == -1: 
          print op
          print sites.cell
          print dot(invcell, newpos), dot(invcell, atom.pos)
          raise internal("Symmetry mapped position not in lattice")
        if index + 1 == label: continue
        intermediate.add(sites[index].label)
      # equivalent atoms share equivalency set.
      for index in intermediate: result[index] = intermediate
  
    return result

  def modisymm(self, structure):
    """ Symmetry modifications, if any. """
    from numpy import allclose
    from .geometry import ModifySymmetry
    if len(self.modifications) == 0: return None
    if self.keepsym: return None
    result = []
    equivs = self._symequivs(structure)
    for label, shell in self.modifications.iteritems():
      if any(label in u for u in result): continue
      indices, split = [shell], {0: [label]}
      for item in equivs[label]:
        if item not in self.modifications: 
          if None in split: split[None].append(item)
          else: split[None] = [item]
          continue
        oshell = self.modifications[item]
        found = False
        for i, cmp in enumerate(indices):
          if allclose(cmp, oshell, atol=1e-8): found = True; break
        if not found: 
          indices.append(shell)
          split[len(indices)-1] = [item]
        elif item not in split[i]: split[i].append(item)
      result.extend(split.values())
    return ModifySymmetry(*result)

  def __get__(self, owner, instance=None): return self
  def __set__(self, owner, value):
    if value is None: self.modifications = {}; return
    self.modifications = value
  def output_map(self, crystal=None, structure=None, **kwargs):
    if len(self.modifications) == 0: return None
    from operator import itemgetter
    from numpy import allclose
    from ..tools.input import Tree
    equivs = self._symequivs(structure)
    if self.keepsym:
      shells = self.modifications.copy()
      # Adds equivalent shells 
      for label, shell in self.modifications.iteritems():
        for equiv in equivs[label]:
          if equiv == label: continue
          elif equiv in shells:
            if not allclose(shell, shells[equiv], atol=1e-8):
              raise ValueError( 'Found equivalent atoms ({0}, {1}) '           \
                                'with different input shells.'                 \
                                .format(label, equiv) )
          else: shells[equiv] = shell
            
          shells[equiv] = shell
    else: shells = self.modifications
    raw = "{0}\n".format(len(shells))
    for key, value in sorted(shells.items(), key=itemgetter(0)):
      raw += "{0}  {1}\n".format(key, " ".join(str(u) for u in value))
    if kwargs.get('breaksym', True) != self.breaksym:
      results = Tree() 
      results['keepsymm' if self.keepsym else 'breaksym'] = True
      results[self.keyword] = raw
      return results
    return {self.keyword: raw}

  def read_input(self, tree, owner=None, **kwargs):
    from numpy import array
    from ..error import ValueError
    if not hasattr(tree, 'prefix'): raise ValueError('Unexpected empty node.')
    self.breaksym = kwargs.get('breaksym', True)

    lines = tree.value.splitlines()
    N = int(lines[0].rstrip().lstrip())
    # breaksym is always false if reading from input.
    self.breaksym = False
    self.modifications = {}
    for line in self.lines[1:N+1]:
      line = line.split()
      self.modifications[int(line[0])] = array(line[1:], dtype='float64')

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Dumps object to string, prettily. """
    result = {}
    if defaults is None:
      result['{0}.breaksym'.format(name)] = repr(self.breaksym)
    elif self.breaksym == False: 
      result['{0}.breaksym'.format(name)] = repr(True)
    for key, value in self.modifications.iteritems():
      result['{0}[{1}]'.format(name,key)] = repr(value) 
    return result

class BasisSet(AttrBlock):
  """ Basis set block.
  
  
      Sub-block defining the basis set of the CRYSTAL_ calculation. In
      practice, it is composed of a dictionary which holds the description of
      each specie-dependent sets of guassians, and of keyword parameters, such
      as :py:attr:`ghosts`.

      Currently, only the general basis sets have been implemented, via the
      :py:class:`Shell`. The input looks as follows:

      .. code-block:: python

        functional.basis['Si'] = [ Shell('s', a0=[16120.0,  0.001959],
                                              a1=[ 2426.0,  0.01493], 
                                              a2=[  553.9,  0.07285],
                                              a3=[  156.3,  0.2461], 
                                              a4=[   50.07, 0.4859],
                                              a5=[   17.02, 0.325]),
                                   Shell('sp', a0=[292.7,   -0.002781, 0.004438],
                                               a1=[ 69.87,  -0.03571,  0.03267],
                                               a2=[ 22.34,  -0.115,    0.1347],
                                               a3=[  8.15,   0.09356,  0.3287],
                                               a4=[  3.135,  0.603,    0.4496]), 
                                   Shell('sp', 4.0, a0=[1.22, 1.0, 1.0]),
                                   Shell('sp', 0.0, a0=[0.55, 1.0, 1.0]),
                                   Shell('sp', 0.0, a0=[0.27, 1.0, 1.0]) ]

         functional.basis.atomsymm = True

      There are few moving parts to the code above.
      
      *Indexing*: The basis set of each atomic function can be accessed
      *equivalently* via its atomic symbol ('Si'), its name ('Silicon', first
      letter capitalized), or its atomic number. To create basis-sets with
      numbers greater than 100, then one should use the integer format only
      (e.g.  201, 301, ...) In practice, the indexing codes for same
      information as the first record of a basis set description.

      *List  of orbitals*: Each item of the list is a single set of orbitals.
      At present, only the general basis set has been implemented. The full
      description is given in :py:class:`Shell`.

      *Usage and keywords*:Unlike CRYSTAL_, any number of species can be
      defined. Which should actually enter the input is determined at run-time.
      This makes it somewhat easier to define functionals for a whole set of
      calculations.  Finally, other parameters from the basis-block are defined
      simply as attributes (see last line of code above). As with other
      input-blocks, additional keywords can be defined via
      :py:meth:`add_keyword`.
  """
  __ui_name__ = 'basis'
  """ Name used when printing user-friendly repr. """
  def __init__(self):
    """ Creates basis set block. """
    from .input import BoolKeyword
    super(BasisSet, self).__init__()
    self._functions = {}
    """ Dictionary holding basis functions. """
    self.atomsymm = BoolKeyword()
    """ Prints point symmetry for each atomic position. """
    self.symmops  = BoolKeyword()
    """ Prints point symmetry operators/ """
    self.charged  = BoolKeyword()
    """ Allows charged system. """
    self.noprint  = BoolKeyword()
    """ Disable basis set printing. """
    self.paramprt = BoolKeyword()
    """ Print code dimensions parameters. """
    self.ghosts   = Ghosts()
    """ Removes atoms but keeps basis functions. 
    
        This keyword should be in the basis subsection of the input.
        In practice, it can be set as follows:
    
        >>> functional.basis.ghosts = [1, 3, 5], False
    
        The second item means that symmetries are *kept*.
        This results in the following crystal input:
    
          | KEEPSYMM
          | GHOSTS
          | 3
          | 1 3 5
    
        The first item in the tuple is a sequence of atomic lables (e.g.
        integers). The second is whether or not to break symmetries. It is True
        by default. In the above, symmetries are kept. Hence atoms equivalent
        those explicitely given are also removed.
    
        The two attributes can also be obtained and set using
        :py:attr:`~Ghosts.value` and :py:attr:`~Ghosts.breaksym`:
    
        >>> functional.basis.ghosts.breaksym
        True
        >>> functional.basis.ghosts.value
        [1, 3, 5], 
    
    """
    self.chemod = Chemod()
    """ Sets initial electronic occupation. """

  @property
  def raw(self):
    """ Raw basis set input. """
    from .. import periodic_table as pt
    # figure out set of species
    # Now print them
    result = ''
    for key, value in self.iteritems():
      if len(value) == 0: continue
      if isinstance(key, str): key = getattr(pt, key).atomic_number
      result += '{0} {1}\n'.format(key, len(value))
      for function in value:
        result += function.raw.rstrip()
        if result[-1] != '\n': result += '\n'
      result = result.rstrip()
      if result[-1] != '\n': result += '\n'
    result += '99 0\n'
    return result
  @raw.setter
  def raw(self, value):
    """ Sets basis data from raw CRYSTAL input. """
    from ..error import NotImplementedError as NA, IOError
    # clean up all basis functions. 
    for key in self.keys(): del self[key]
    lines = value.split('\n')
    while len(lines):
      line = lines.pop(0).split()
      type, nshells = int(line[0]), int(line[1])
      if type == 99: break
      type = specie_name(type)
      self._functions[type] = []
      for i in xrange(nshells):
        bt = int(lines[0].split()[0])
        if bt == 0:
          ncoefs = int(lines[0].split()[2])
          if len(lines) < ncoefs + 1: raise IOError('Unexpected end-of-file.')
          basis = Shell()
          basis.raw = '\n'.join(lines[:ncoefs+1])
          self._functions[type].append(basis)
          lines = lines[ncoefs+1:]
        else:
          raise NA('Basis type {0} hasnot yet been implemented.'.format(bt))


  def __setitem__(self, specie, shell):
    """ Sets basis function for a particular specie. """
    specie = specie_name(specie)
    self._functions[specie] = shell
  def __getitem__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = specie_name(specie)
    return self._functions[specie]
  def __delitem__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = specie_name(specie)
    del self._functions[specie]
  def __contains__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = specie_name(specie)
    return specie in self._functions
  def __len__(self):
    """ Number of atomic species. """
    return len(self._functions)
  def __iter__(self):
    """ Iterates over species. """
    for key in self._functions: yield key
  def keys(self):
    """ List of species. """
    return self._functions.keys()
  def items(self):
    """ Tuples (specie, basis-functions) """
    return self._functions.items()
  def values(self):
    """ List of basis functions, """
    return self._functions.values()
  def iterkeys(self):
    """ Iterates over species. """
    return self._functions.iterkeys()
  def iteritems(self):
    """ Iterates over tuple (key, basis-functions). """
    return self._functions.iteritems()
  def itervalues(self):
    """ Iterates over basis-functions. """
    return self._functions.itervalues()

  def append(self, specie, shell):
    """ Adds shell to set. """
    specie = specie_name(specie)
    # now adds shell correctly.
    if specie in self._functions: self._functions[specie].append(shell)
    else: self._functions[specie] = [shell]
    return self

  def __repr__(self, defaults=False, name=None):
    """ Returns representation of this instance """
    from ..tools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def output_map(self, **kwargs):
    """ Dumps CRYSTAL input to string. """
    if kwargs.get('structure', None) is not None:
      from copy import deepcopy
      c = deepcopy(self)
      species = set([a.type for a in kwargs['structure'].eval()])
      for specie in species:
        if specie not in c:
          raise KeyError('No basis set for specie {0}.'.format(specie))
      for specie in self.keys():
        if specie not in species: del c[specie]
    else: c  = self
    results = super(BasisSet, self).output_map(**kwargs)
    results.prefix = c.raw
    return results

  def read_input(self, tree, owner=None, **kwargs):
    """ Parses an input tree. """
    self._functions = {}
    self._input = self.__class__()._input
    super(BasisSet, self).read_input(tree, owner=self, **kwargs)
    if hasattr(tree, 'prefix'): self.raw = tree.prefix

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Prettier representation """
    from pprint import pprint
    from StringIO import StringIO
    from ..tools.uirepr import add_to_imports
    results = super(BasisSet, self).__ui_repr__( imports, name, defaults,
                                                 ['raw', 'gaussian94'] )
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    for key, value in self._functions.iteritems():
      newname = '{0}[{1!r}]'.format(name, key)
      string = StringIO()
      pprint(value, string)
      string.seek(0)
      results[newname] = string.read()
      for u in value: add_to_imports(u, imports)
    return results
  
  def update(self, basis):
    """ Updates dictionary of basis-sets.
    
        Overwrites existing species.
    """
    for key, value in basis.iteritems():
      if key in self: del self[key]
      for v in value: self.append(key, v)

  @property
  def gaussian94(self):
    """ Basis set in GAUSSIAN94 format """
    result = '\n\n****\n'
    for specie, basis in self._functions.iteritems():
      result += '{0} 0\n'.format(specie)
      for shell in basis: 
        result += '{0} {1} 1.0\n'.format(shell.type.upper(), len(shell))
        for function in shell.functions:
          alpha, coef1 = function[0], function[1]
          if hasattr(alpha, 'rescale'):
            alpha = float(alpha.rescale(crystal_invbohr2).magnitude)
          if len(function) == 2:
            result += '{0:> 20.12f}    {1:> 20.12f}\n'.format(alpha, coef1)
          else:
            result += '{0:> 20.12f}    {1:> 20.12f}  {2[2]:> 20.12f}\n'               \
                      .format(alpha, coef1, function)
      result += '****\n'
    return result

  def add_specie(self, string):
    """ Add specie using CRYSTAL string input.
    
        Just a convenience function to quickly add species by copy/pasting
        crystal input. 

        :param str string: 
           Full CRYSTAL_ input for an atomic specie.
    """
    a = BasisSet()
    try: a.raw = string
    except: 
      from sys import exc_info
      type, value, traceback = exc_info()
      message = string + '\n\n ******* Incorrect input.' 
      if value is not None: value = value, string
      else: type.args = tuple(list(type.args) + [message])
      raise type, value, traceback
    self.update(a)
     

  

