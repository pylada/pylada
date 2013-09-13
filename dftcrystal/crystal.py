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
__all__ = ['Crystal']
from .molecule import Molecule
class Crystal(Molecule):
  """ CRYSTAL-wise structure, e.g. functional approach. 

      CRYSTAL_ proposes a functional_ approach to crystals, as opposed to the
      imperative_ style used by Pylada. In practice, this means that CRYSTAL_
      declares a chain of functions which acts upon an initial data and
      transform it. The data is generally a space-group with a set of atomic
      sites. The functions can be affine transformations on these sites,
      additions and removals of sites, transformation upon the space-group, and
      even strain relaxation.

      In practice, both approaches can "do the same thing", and indeed both
      CRYSTAL_ and Pylada provide similar functionalities, e.g. creating a
      supercell from  a unit-cell. However, there are clear benefits from
      allowing users to keep working the CRYSTAL_ way when working with
      CRYSTAL_. 

      This class provides a python wrapper around CRYSTAL_'s approach. It is
      only a wrapper, in that it merely contains the data necessary to create a
      CRYSTAL_ input file. It does not re-implement any of CRYSTAL_'s
      algorithms. However, it does provides the ability to call CRYSTAL_ and
      extract the results as a Pylada
      :py:class:`~pylada.crystal.cppwrappers.Structure` instance. In this way,
      both approaches can be mixed, allowing for a complete integration of
      CRYSTAL_ with python. At present, both the initial data and chain of
      functions are represented and stored within the same :py:class:`Crystal`
      instance. This is a bit sub-optimal. A better implementation would make
      clear the separation bewteen initial data and functions. However, it is
      less wordy.

      The initial data upon which CRYSTAL_ applies its chain of functions is
      defined using the :py:attr:`symmgroup`, :py:attr:`args`, and
      :py:attr:`atoms`. The first references the space-group of the crystal,
      the second defines the lattice parameters, and the last holds the
      (initial) atomic sites.  :py:class:`Crystal` is also a list where each
      item is an operation upon the initial crystal. 

      In practice, the initial data is declared as follows:

      >>> from pylada.dftcrystal import Crystal
      >>> crystal = Crystal(227, 5.43)                                        \\
      ...                  .add_atom(0.125, 0.125, 0.125, "Si")

      This declares the diamond silicon structure. The first argument is the
      space-group, and the second the only required lattice parameter. If other
      parameters were needed they would be listed directly ``Crystal(?, a, b,
      c, alpha, beta, gamma)``. Only the parameters required for that
      particular space-group should be listed.

      The second line adds an atomic site to the initial structure.
      :py:meth:`add_atom` returns an the :py:class:`Crystal` instance itself,
      so that atomic site declarations can be chained:

      >>> crystal = Crystal(...)                                              \\
      ...           .add_atom(...)                                            \\
      ...           .add_atom(...)                                            \\
      ...           .add_atom(...)
      
      Only those atoms which are strictly inequivalent by symmetry operations
      should be listed.

      :py:class:`Crystal` instances function as lists of transformations which
      are to be applied to the initial structure. The simplest way to add an
      operation is to use the traditional CRYSTAL_ format:


      >>> crystal = Crystal(227, 5.43)                                        \\
      ...                  .add_atom(0.125, 0.125, 0.125, "Si")               \\
      ...                  .append('supercel', '2 0 0 0 1 0 0 0 1')

      Here we have created a supercell of diamond with two unit cell in the
      (100) direction. The first string in :py:meth:`append` is the keyword and
      the second string its input. Much as :py:meth:`add_atom`,
      :py:meth:`append` returns the calling instance of :py:class:`Crystal` so
      that calls can be chained into a single declarative expression.

      A number of operations are implemented in a more pythonic manner. These
      can be added to the chain of functions by calling :py:meth:`append` with
      an operation instance as a the only argument.

      >>> from pylada.dftcrystal import Slabcut
      >>> crystal.append( Slabcut(hkl=(1, 0, 0), isup=1, nl=3) )

      :py:class:`~pylada.dftcrystal.input.Slabcut` is an operation to create a
      thin-film from a 3d bulk material. 

      Finally, the whole  "data+functions" object can be evaluated with
      :py:meth:`eval`. This will return a
      :py:class:`~pylada.crystal.cppwrappers.Structure` instance which can be
      used with other Pylada functionalities. Internally, :py:meth:`eval` makes a
      call to CRYSTAL_ and greps the output to construct the output structure.

      :param symmgroup:
         String or integer denoting the space-group.
      :param args:
         Lattice parameters :math:`a`, :math:`b`, :math:`c`, :math:`\\alpha`,
         :math:`\\beta`, :math:`\\gamma`. In practice, these parameters are
         listed as given in the third (or fourth) line of the input to
         CRYSTAL_. Hence it should conform to the same syntax and the same
         order.
      :param shift: 
         Optional. Non-standard shift of the origin. Should be a sequence of
         three floats.

      .. seealso::
      
         For operations, see :py:mod:`~pylada.dftcrystal.geometry`.


      .. _functional: http://en.wikipedia.org/wiki/Functional_programming
      .. _imperative: http://en.wikipedia.org/wiki/Imperative_programming
  
  """
  keyword = 'crystal'
  """ CRYSTAL keyword. """
  def __init__(self, symmgroup=1, *args, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    shift = kwargs.pop('shift', 0)
    super(Crystal, self).__init__(symmgroup, **kwargs)

    self.params = list(args)
    """ Lattice parameters. """
    self.shift = shift
    """ Non-standard shift of the origin. """

  @property
  def shift(self):
    """ Shift of the origin of the cell,
    
         - if 0 or second, second setting given international tables
           (equivalent to IFSO == 0).
         - if 1 or first, second setting given international tables
           (equivalent to IFSO == 1).
         - if a sequence of three, then an actual shift in reciprocal space * 24.
    """
    return self._shift
  @shift.setter
  def shift(self, value):
    from numpy import array
    from ..error import ValueError
    if hasattr(value, '__iter__'): 
      self._shift = array(value, dtype='float64')
    elif value == 0 or value == 'second': self._shift = 0
    elif value == 1 or value == 'first': self._shift = 1
    else: raise ValueError('Unknown input value for shift.')

  def __repr_header__(self):
    """ Dumps representation to string. """
    # prints construction part.
    length = len('{0.__class__.__name__}('.format(self)) + 1
    args = [repr(self.symmgroup)]
    for o in self.params: args.append('{0!r}'.format(o))
    indent = ' '.join('' for i in xrange(length))
    for key, value in self.__dict__.iteritems():
      if key in ['atoms', 'symmgroup', 'params'] : continue
      if key == '_shift': key = 'shift'
      args.append('\\\n{2}{0}={1!r}'.format(key, value, indent))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  @property 
  def raw(self):
    """ Raw input for structure. """
    from ..physics.spacegroup import space_group as sg
    from ..error import ValueError
    from ..periodic_table import find as find_specie

    if hasattr(self.shift, '__iter__'): hasshift = 2
    else:                               hasshift = self.shift

    try: symmgroup = int(self.symmgroup)
    except: 
      if self.symmgroup not in sg:
        raise ValueError('Unknown space group {0.symmgroup}'.format(self))
      symmgroup = sg[self.symmgroup].index
    
    # control flags + symmgroup 
    result = '0 {0} {1}\n{2}\n'.format( getattr(self, 'ifhr', 0),              \
                                        hasshift, symmgroup )
    # shift if needed, watch out for weird ass input.
    if hasshift == 2:
      shift = self.shift * 24.0 + 0.1
      shift = [int(u) for u in shift]
      result += '{0[0]} {0[1]} {0[2]}\n'.format(shift)

    # lattice parameters.
    for o in self.params: result += '{0} '.format(o)

    # number of atoms + atoms
    result += '\n{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
      try: n = find_specie(name=atom.type).atomic_number
      except:
        try: n = int(atom.type)
        except: 
          raise ValueError( 'Could not make sense of atomic type {0.type}.'    \
                            .format(atom) )
      result += '{0} {1.pos[0]} {1.pos[1]} {1.pos[2]}\n'.format(n, atom)
    return result

  @raw.setter
  def raw(self, value):
    """ Reads crystal input. """
    from ..periodic_table import find as find_specie
    from numpy import array
    if not hasattr(value, '__iter__'): value = value.split('\n')
    data = value.pop(0).split()
    self.symmgroup, self.ifhr, shift = (int(u) for u in data[:3]) 
    data = value.pop(0)
    if self.symmgroup == 0: self.symmgroup = int(data.split()[0])
    else: self.symmgroup = data

    if shift > 1: 
      self.shift = array(value.pop(0).split(), dtype='float64') / 24.0
    else: self.shift = shift

    self.params = [float(u) for u in value.pop(0).split()]
    n = int(value.pop(0).split()[0])
    self.atoms = []
    for line in value[:n]:
      line = line.split()
      type = int(line[0])
      if type < 100: type = find_specie(atomic_number=type).symbol
      self.add_atom(pos=[float(u) for u in line[1:4]], type=type)
