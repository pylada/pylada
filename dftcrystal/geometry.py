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
__all__ = [ 'RemoveAtoms', 'Slabinfo', 'Slabcut', 'Slab',
            'DisplaceAtoms', 'InsertAtoms', 'ModifySymmetry',
            'AffineTransform', 'Elastic', 'Supercell', 'Supercon',
            'KeepSymm', 'BreakSym', 'BohrUnits', 'AngstromUnits',
            'Fractional' ]
from ..tools.input import BaseKeyword, BoolKeyword

class KeepSymm(BoolKeyword):
  """ Switches to keeping symmetries. """
  keyword = "keepsymm"
  def output_map(self, **kwargs):
    """ Prints to Tape.

        Does not print if breaksym is False on input.
    """
    if kwargs.get('breaksym', True) == False: return None
    return {self.keyword: True}
class BreakSym(BoolKeyword):
  """ Switches to breaking symmetries (Default). """
  keyword = "breaksym"
  def output_map(self, **kwargs):
    """ Prints to Tape.

        Does not print if breaksym is True on input.
    """
    if kwargs.get('breaksym', False) == True: return None
    return {self.keyword: True}
class BohrUnits(BoolKeyword):
  """ Switches to Bohr units. """
  keyword = "bohr"
  def output_map(self, **kwargs):
    """ Prints to Tape.

        Does not print if units is bohr on input.
    """
    if kwargs.get('units', '') == 'bohr': return None
    return {self.keyword: True}
class AngstromUnits(BoolKeyword):
  """ Switches to angstrom units (Default). """
  keyword = "angstrom"
  def output_map(self, **kwargs):
    """ Prints to Tape.

        Does not print if units is angstrom on input.
    """
    if kwargs.get('units', '') == 'angstrom': return None
    return {self.keyword: True}
class Fractional(BoolKeyword):
  """ Switches to fractional units. """
  keyword = "fraction"
  def output_map(self, **kwargs):
    """ Prints to Tape.

        Does not print if units is fraction on input.
    """
    if kwargs.get('units', '') == 'fraction': return None
    return {self.keyword: True}


class RemoveAtoms(BaseKeyword):
  """ Remove atoms from structure.
  
      Atoms are distinguished by their CRYSTAL label.
      To remove atoms 1 and 5:

      >>> crystal.append( RemoveAtoms(1, 5) )

      To modify/check that list

      >>> crystal[-1].labels
      [1, 5]
      >>> crystal[-1].labels[1] = 6
      >>> crystal[-1].labels
      [1, 6]
  """
  keyword = 'atomremo'
  """ CRYSTAL keyword. """
  def __init__(self, *args, **kwargs):
    """ Creates atom removal operator.
    
        :param args:
          Atoms to remove are set by specifying them as arguments to the
          iterator.
    """ 
    from copy import copy
    super(RemoveAtoms, self).__init__(**kwargs)
    self.labels = copy(args)
    """ Indices of atoms to remove. """

  @property
  def raw(self):
    """ Returns input string. """
    return str(len(self.labels)) + '\n'                                        \
           + ' '.join(str(u) for u in self.labels) + '\n'
  @raw.setter
  def raw(self, value):
    """ Read crystal input. """
    value = value.split('\n')
    n = int(value[0].split()[0])
    i = 1
    self.labels = []
    while len(self.labels) < n:
      self.labels.extend([u for u in value[i].split()])
      i += 1
    self.labels = [int(u) for u in self.labels[:n]]

  def __repr__(self):
    result = '{0.__class__.__name__}({1}, '                                    \
             .format( self, ', '.join(str(u) for u in self.labels))
    return result[:-2] + ')'
    
class Slabinfo(BaseKeyword):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl=None):
    from numpy import array
    super(Slabinfo, self).__init__()
    self.hkl = hkl
    """ Surface normal. """
    if hkl is not None: self.hkl = array(hkl, dtype='int32')
  @property
  def raw(self):
    """ Input to SLABINFO """
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}'.format(self)
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    self.hkl = array(value.split()[:3], dtype='int32')
  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}])'     \
           .format(self)
  def output_map(self, **kwargs):
    """ Prints input. """
    if self.hkl is None: return None
    return super(Slabinfo, self).output_map(**kwargs)

class Slabcut(Slabinfo):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl=None, isup=None, nl=None, raw=None):
    super(Slabcut, self).__init__(hkl)
    self.isup = isup
    """ Subsurface atomic label. """
    self.nl = nl
    """ Number of monolayers. """
  @property
  def raw(self):
    """ Input to SLABINFO """
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}\n{0.isup} {0.nl}\n'.format(self)
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    line0, line1 = value.split('\n')[:2]
    self.hkl = array(line0.split()[:3], dtype='int32')
    self.isup, self.nl = [int(u) for u in line1.split()[:2]]

  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}], '    \
           '{0.isup}, {0.nl})'.format(self)
Slab = Slabcut
""" Alias for :py:class`~pylada.dftcrystal.geometry.Slabcut`. """

class DisplaceAtoms(BaseKeyword):
  """ Displaces atoms.
  
      This keywords applies a displacement to a set of atoms, identified by
      their labels:

      >>> from pylada.dftcrystal import DisplaceAtoms
      >>> disp = DisplaceAtoms()                                               \\
      ...                     .add_atom(0, 0.01, 0.002, 1)                     \\
      ...                     .add_atom(0, -0.01, 0.05, 5)
      >>> structure.append(disp)


      The above creates a displacement operations for two atoms, labelled 1 and
      5. The displacements are the current units of the crystal structure. It
      would look as follows in CRYSTAL_'s input:

        | KEEPSYMM
        | ATOMDISP
        | 2
        | 1 0.0 0.01 0.002
        | 5 0.0 -0.01 0.05
  """
  keyword = 'atomdisp'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    """ Creates a displacement field. """
    super(DisplaceAtoms, self).__init__(**kwargs)
    self.atoms = []
    """ Atomic displacements. """

  def add_atom(self, *args, **kwargs):
    """ Adds a displacement to a given atom.
    
        At present, atom.type should be an index to an atom in the structure.
    """
    from ..crystal import Atom
    self.atoms.append(Atom(*args, **kwargs))
    return self

  def __repr__(self):
    """ Dumps representation to string. """
    result = super(DisplaceAtoms, self).__repr__()
    indent = ' '.join('' for i in xrange(result.find('(')))
    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent, dummy)
    return result

  @property
  def raw(self):
    """ Raw input to CRYSTAL. """
    result = '{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
      result += '{0.type} {0.pos[0]} {0.pos[1]} {0.pos[2]}\n'.format(atom)
    return result
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    self.atoms = []
    value = value.split('\n')
    n = int(value.pop(0).split()[0])
    for line in value[:n]: 
      line = line.split()
      type = int(line[0])
      pos = array(line[1:4], dtype='float64')
      self.add_atom(type=type, pos=pos)
    
   
class InsertAtoms(DisplaceAtoms):
  """ insert atom into structure. """
  keyword = 'atominse'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    super(InsertAtoms, self).__init__(**kwargs)

  @property
  def raw(self):
    """ Raw input to CRYSTAL. """
    from .. import periodic_table as pt
    # number of atoms + atoms
    result = '{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
      n = getattr(pt, atom.type, None)
      if n is not None: n = n.atomic_number
      else:
        try: n = int(atom.type)
        except: 
          raise ValueError( 'Could not make sense of atomic type {0.type}.'    \
                            .format(atom) )
      result += '{0} {1.pos[0]} {1.pos[1]} {1.pos[2]}\n'.format(n, atom)
    return result
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    from .. import periodic_table as pt
    self.atoms = []
    value = value.split('\n')
    n = int(value.pop(0).split()[0])
    for line in value[:n]: 
      line = line.split()
      type = int(line[0])
      if type < 100: 
        for key, value in pt.__dict__.iteritems():
          if getattr(pt.__dict__[key], 'atomic_number', -1) == type:
            type = key
            break
      pos = array(line[1:4], dtype='float64')
      self.add_atom(type=type, pos=pos)

class ModifySymmetry(BaseKeyword):
  """ Modify symmetries. """
  keyword = 'modisymm'
  """ CRYSTAL keyword. """
  def __init__(self, *args):
    """ Creates symmetry modifier.

        Each argument is a sequence of atomic label. Each argument will be
        assigned a different flag. This input is somewhat more condensed than
        the original CRYSTAL input.

        >>> a = ModifySymmetry([1, 2, 3], [4, 5, 6])
        >>> structure.append(a)
        
	The snippet above adds a symmetry modifier to a structure which
	disables any symmetry operation  between the group of atoms 
	``[1, 2, 3]`` and the group of atoms ``[4, 5, 6]``.
    """
    super(ModifySymmetry, self).__init__()
    self.groups = []
    """ Atoms for which to modify symmetries. """
    for o in args:
      self.groups.append(list(o) if hasattr(o, '__iter__') else [o])

  @property
  def raw(self):
    """ Raw CRYSTAL input. """
    result = '{0}\n'.format(sum(len(o) for o in self.groups))
    for i, labels in enumerate(self.groups):
      result += ' '.join('{0} {1}'.format(u, i+1) for u in labels) + '\n'
    return result
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    value = value.split()
    n = int(value.pop(0))
    d = {}
    for label, flag in zip(value[:2*n:2], value[1:2*n+1:2]):
      if flag in d: d[flag].append(int(label))
      else: d[flag] = [int(label)]
    self.groups = d.values()
  def __repr__(self):
    """ Representation of this instance. """
    return '{0.__class__.__name__}({1})'                                       \
           .format(self, ', '.join(repr(u) for u in self.groups))

class AffineTransform(BaseKeyword):
  """ Affine transformation applied to the crystal or crystal fragment. """
  keyword = 'atomrot'
  """ CRYSTAL keyword. """
  def __init__( self, labels=None, vectrans=None, origin=None, bondtrans=None,
                euler=None, bondrot=None, bondtoz=None):
    """ Creates rotation. """
    self.labels = labels
    """ Selects what to apply transform to. 

        Should be either:

          - 'all' or None or an empty sequence, for the whole crystal.
          - a sequence of integers signifying the atoms to rotate
          - a single integer or a sequence of one integer, signifying that the
            molecule to which this atom belongs will be rotated.
    """
    from ..error import ValueError
    super(AffineTransform, self).__init__()

    if [vectrans is None, origin is None, bondtrans is None].count(False) > 1:
      raise ValueError('More than one type of translation was selected.')
    if [euler is None, bondrot is None, bondtoz is None].count(False) > 1:
      raise ValueError('More than one type of rotation was selected.')

    self.vectrans  = vectrans
    self.bondtrans = bondtrans
    self.origin    = origin
    self.euler     = euler
    self.bondrot   = bondrot
    self.bondtoz   = bondtoz

  @property
  def vectrans(self):
    """ Selects translation by a vector. 

        Should be None (unselected) or a sequence of three real numbers
        defining the translation.
    """ 
    return self._vectrans
  @vectrans.setter
  def vectrans(self, value): 
    from numpy import array
    from quantities import angstrom
    if value is None: self._vectrans = None; return
    if hasattr(value, 'rescale'): self._vectrans = value 
    else: self._vectrans = array(value, dtype='float64') * angstrom
    self._origin, self._bondtrans = None, None
  @property
  def bondtrans(self):
    """ Selects translation by a vector. 

        Should be None (unselected) or a sequence of two labels(integers) and a
        real number. The first two indicate the axis and direction of the
        translation (first -> second), and the third its modulus in angstrom.

        .. _quantities: http://packages.python.org/quantities/index.html
    """ 
    return self._bondtrans
  @bondtrans.setter
  def bondtrans(self, value): 
    self._bondtrans = value
    if value is not None:
      self._origin, self._vectrans = None, None
  @property
  def origin(self):
    """ Shifts origin to specific atom. 
    
        Can be None or an atom label.
    """
    return self._origin
  @origin.setter
  def origin(self, value):
    if value is None: self._origin = None; return
    self._origin = int(value) 
    self._bondtrans, self._vectrans = None, None

  @property
  def euler(self):
    """ Defines a rotation using Euler matrices and a given atomic origin. 
    
        Should consist of three real numbers defining the Euler rotations (in
        degrees), and an atomic index defining the origin of the cartesian
        system. 
    """
    return self._euler
  @euler.setter
  def euler(self, value):
    self._euler = value
    if value is not None:
      self._bondrot, self._bondtoz = None, None
  @property
  def bondrot(self):
    """ Defines a rotation via two atoms forming the axis. 
 
        Should be two atomic indices followed by the rotation angle in degrees.
    """
    return self._bondrot
  @bondrot.setter
  def bondrot(self, value):
    self._bondrot = value
    if value is not None:
      self._euler, self._bondtoz = None, None
  @property
  def bondtoz(self):
    """ Defines a rotation via two atoms forming the axis. 
 
        Should be two atomic indices followed by the rotation angle in degrees.
    """
    return self._bondtoz
  @bondtoz.setter
  def bondtoz(self, value):
    self._bondtoz = value
    if value is not None:
      self._euler, self._bondrot = None, None

  @property
  def raw(self):
    """ Creates raw input for crystal. """
    from quantities import angstrom, degree
    # first line
    if self.labels is None: result = '0\n'
    elif not hasattr(self.labels, '__iter__'): result = str(-self.labels) + '\n'
    elif len(self.labels) == 0: result = '0\n'
    elif len(self.labels) == 1: result = str(-self.labels[0]) + '\n'
    else:
      result = '{0}\n{1}\n'.format( len(self.labels),                          \
                                    ' '.join(str(u) for u in self.labels) )
    if self.origin is not None: result += str(self.origin) + ' '
    elif self.bondtrans is not None: result += '0 '
    elif self.vectrans is not None: result += '-1 '
    else: result += '999 '

    if self.bondtoz is not None: result += '1\n'
    elif self.bondrot is not None: result += '1\n'
    elif self.euler is not None: result += '-1\n'
    else: result += '999\n'

    if self.bondtrans is not None:
      a, b, mod = int(self.bondtrans[0]), int(self.bondtrans[1]), self.bondtrans[2]
      if hasattr(mod, 'rescale'): mod = float(mod.rescale(angstrom).magnitude)
      result += '{0} {1} {2}\n'.format(a, b, mod)
    elif self.vectrans is not None:
      vec = self.vectrans
      if hasattr(vec, 'rescale'): vec = vec.rescale(angstrom).magnitude
      result += '{0[0]} {0[1]} {0[2]}\n'.format(vec)

    if self.bondtoz is not None: 
      a, b = int(self.bondtoz[0]), int(self.bondtoz[1])
      result += '{0} {1} 0\n'.format(a, b)
    elif self.bondrot is not None:
      a, b = int(self.bondrot[0]), int(self.bondrot[1])
      c = self.bondrot[2]
      if hasattr(c, 'rescale'): c = int(c.rescale(degree).magnitude + 0.01)
      result += '{0} {1} {2}\n'.format(a, b, c)
    elif self.euler is not None: 
      a, b, c = self.euler[:3]
      if hasattr(a, 'rescale'): a = int(a.rescale(degree).magnitude + 0.01)
      if hasattr(b, 'rescale'): b = int(b.rescale(degree).magnitude + 0.01)
      if hasattr(c, 'rescale'): c = int(c.rescale(degree).magnitude + 0.01)
      d = int(self.euler[3])
      result += '{0} {1} {2} {3}\n'.format(a, b, c, d)
    return result

  @raw.setter
  def raw(self, value):
    """ Read crystal input. """ 
    from numpy import array
    from quantities import angstrom, degree
     
    self.vectrans, self.bondtrans, self.origin = None, None, None
    self.euler, self.bondrot, self.bondtoz = None, None, None

    value = value.split()
    line = int(value.pop(0))
    if line == 0: self.labels = None
    elif line > 0: 
      self.labels = [int(u) for u in value[:line]]
      value = value[line:]
    else: self.labels = line
    trans, rot = [int(u) for u in value[:2]]
    value = value[2:]

    if trans > 0 and trans != 999: self.origin = trans
    elif trans == 0:
      self.bondtrans = int(value[0]), int(value[1]), float(value[2]) * angstrom
      value = value[3:]
    elif trans != 999: 
      self.vectrans = array(value[:3], dtype='float64') * angstrom        
      value = value[3:]
    
    if rot < 0: 
      self.euler = float(value[0]) * degree,                                   \
                   float(value[1]) * degree,                                   \
                   float(value[2]) * degree,                                   \
                   int(value[3])
    elif rot != 999:
      a, b, alpha = int(value[0]), int(value[1]), int(value[2])
      if alpha == 0: self.bondtoz = a, b
      else: self.bondrot = a, b, alpha * degree

  def __repr__(self): 
    """ Representation of AffineTransform. """
    args = []
    if self.labels is None: args.append('labels=None')
    elif isinstance(self.labels, int):
      args.append('labels={0.labels}'.format(self))
    elif len(self.labels) == 0: args.append('labels=None')
    elif len(self.labels) == 1:
      args.append('labels={0.labels[0]}'.format(self))
    else: args.append('labels={0.labels}'.format(self))
    if self.vectrans is not None:
      args.append('vectrans={0!r}'.format(self.vectrans))
    elif self.bondtrans is not None:
      args.append('bondtrans={0!r}'.format(self.bondtrans))
    elif self.origin is not None:
      args.append('origin={0!r}'.format(self.origin))
    if self.euler is not None:
      args.append('euler={0!r}'.format(self.euler))
    elif self.bondrot is not None:
      args.append('bondrot={0!r}'.format(self.bondrot))
    elif self.bondtoz is not None:
      args.append('bondtoz={0!r}'.format(self.bondtoz))
    return "{0.__class__.__name__}(".format(self) + ', '.join(args) + ')'

class Elastic(BaseKeyword):
  """ Elastic deformation of the lattice """
  keyword = 'elastic'
  """ CRYSTAL keyword """
  def __init__(self, matrix=None, is_epsilon=True, constvol=False, **kwargs):
    """ Creates cell-shape deformation. """
    from numpy import array
    super(Elastic, self).__init__(**kwargs)
    self.matrix = matrix
    """ Deformation matrix. 

        It may or may not include the identity, depending on
        :py:attr:`is_epsilon`.
    """
    if matrix is not None: self.matrix = array(matrix)
    self.is_epsilon = is_epsilon
    """ Whether the matrix includes the identity. """
    self.const_volume = constvol
    """ Whether this is a constant volume deformation. """
  @property
  def raw(self):
    if self.matrix is None: return ""
    type = 2 if self.is_epsilon else 1
    if not self.const_volume: type = -type
    result = str(type) + '\n'
    if self.is_epsilon:
      for i in xrange(3):
        result += ' '.join(str(self.matrix[i,j]) for j in xrange(3)) + '\n'
    else: 
      for i in xrange(3):
        result += ' '.join(str(self.matrix.T[i,j]) for j in xrange(3)) + '\n'
    return result
  @raw.setter
  def raw(self, value):
    from numpy import array
    value = value.split('\n')
    type = int(value[0])
    self.is_epsilon = abs(type) == 2
    self.const_volume = type > 0
    self.matrix = array([u.split()[:3] for u in value[1:4]], dtype='float64')
    if not self.is_epsilon: self.matrix = self.matrix.T

  def __repr__(self):
    """ Representation of this object. """
    args = []
    if self.matrix is not None:
      args.append(repr(self.matrix.tolist()))
    if not self.is_epsilon:
      args.append('is_epsilon=False' if len(args) > 0 else 'False')
    if self.const_volume:
      args.append('constvol=True' if len(args) > 0 else 'True')
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))


  def output_map(self, **kwargs):
    if self.matrix is None: return None
    return super(Elastic, self).output_map(**kwargs)

class Supercell(BaseKeyword):
  """ CRYSTAL function to create a supercell. """
  keyword = "supercel"
  def __init__(self, matrix=None):
    from numpy import array
    from ..error import ValueError
    super(Supercell, self).__init__()
    self.matrix = matrix
    if self.matrix is not None:
      self.matrix = array(matrix, dtype='float64')                       \
                            .round().astype('int64')
      N = len(self.matrix.flat)
      if N == 1: self.matrix = self.matrix.reshape(1,1)
      elif N == 4: self.matrix = self.matrix.reshape(2,2)
      elif N == 9: self.matrix = self.matrix.reshape(3,3)
      else: raise ValueError('Incorect matrix size in Supercell.')
    
  @property 
  def matrix(self):
    """ Supercell in basis of the unit-cell. """
    return self._matrix
  @matrix.setter
  def matrix(self, value):
    from numpy import array
    from numpy.linalg import det
    from ..error import ValueError
    if value is None: self._matrix = None; return
    self._matrix = array(value, dtype='float64').round().astype('int64')
    N = len(self._matrix.flat)
    if N == 1: self._matrix = self._matrix.reshape(1,1)
    elif N == 4: self._matrix = self._matrix.reshape(2,2)
    elif N == 9: self._matrix = self._matrix.reshape(3,3)
    else: raise ValueError('Incorect matrix size in Supercell.')
    if det(self._matrix) == 0:
      raise ValueError('Determinant of supercell cannot be zero.')

  def __repr__(self):
    """ Dumps object to stream. """
    if self.matrix is None: return '{0.__class__.__name__}()'.format(self)
    return '{0.__class__.__name__}({1!r})'                                     \
           .format(self, self.matrix.tolist())

  def output_map(self, structure=None, **kwargs):
    """ Checks dimensionality. """
    from numpy import any, abs
    from ..error import ValueError
    if self.matrix is None: return None
    # dimensionality check, wherever possible.
    if structure is not None: 
      for i, op in enumerate(structure): 
        if op is self: break
      if op is self: 
        structure = structure.copy()
        structure[:] = structure[:i]
        try: cell = structure.eval().cell
        except: pass
        else: 
          N = len(self.matrix)
          if N < 1 and any(abs(cell[:, 0] - [5e2, 0,0]) > 1e-8):
            raise ValueError( 'Structure and supercell '                       \
                              'dimensionality do not match.' )
          if N < 2 and any(abs(cell[:, 1] - [0, 5e2, 0]) > 1e-8):
            raise ValueError( 'Structure and supercell '                       \
                              'dimensionality do not match.' )
          if N < 3 and any(abs(cell[:, 2] - [0, 0, 5e2]) > 1e-8):
            raise ValueError( 'Structure and supercell '                       \
                              'dimensionality do not match.' )

    string = '\n'.join(' '.join(str(i) for i in j) for j in self.matrix)
    return {self.keyword: string}
             
  def read_input(self, tree, owner=None, **kwargs):
    from numpy import array
    self.matrix = array(tree.split(), dtype='float64')                         \
                       .round().astype('int64')

class Supercon(Supercell):
  keyword = 'supercon'
