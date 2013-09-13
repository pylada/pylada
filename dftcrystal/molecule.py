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
__all__ = ['Molecule']
from ..tools.input import ListBlock, BaseKeyword
class Molecule(ListBlock):
  """ Base class for functional crystal structures for CRYSTAL. """
  keyword = 'molecule'
  def __init__(self, symmgroup=1, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    super(Molecule, self).__init__()

    self.symmgroup = symmgroup
    """ Symmetry group. """
    self.atoms = []
    """ List of atoms. """
    for key, value in kwargs.iteritems(): setattr(self, key, value)

  def add_atom(self, *args, **kwargs):
    """ Adds an atom to the structure. """
    from ..crystal import Atom
    atom = Atom(*args, **kwargs)
    self.atoms.append(atom)
    return self

  def __repr_header__(self):
    """ Dumps representation header to string. 
    
        Mostly useful to limit how much should be rewritten of repr in derived
        classes.
    """
    # prints construction part.
    length = len('{0.__class__.__name__}('.format(self)) + 1
    args = [repr(self.symmgroup)]
    indent = ''.join([' ']*length) 
    for key, value in self.__dict__.iteritems():
      if key in ['atoms', 'symmgroup'] : continue
      args.append('\\\n{2}{0}={1!r}'.format(key, value, indent))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Dumps representation to string. """
    from ..tools.uirepr import add_to_imports
    result = self.__repr_header__()
    indent = ' '.join('' for i in xrange(result.find('(')+1))
    add_to_imports(self, imports)

    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent, dummy)

    for item in self: 
      if item.__class__ is BaseKeyword:
        args = [repr(item.keyword)]
        if getattr(item, 'raw', None) is not None: args.append(repr(item.raw))
        result += '\\\n{0}.append({1})'.format(indent, ', '.join(args)) 
      else:
        add_to_imports(item, imports)
        item = repr(item).rstrip().lstrip()
        item = item.replace('\n', '\n' + indent)
        result += '\\\n{0}.append({1})'.format(indent, item) 

    result = result.rstrip()
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    return {name: result.rstrip()}

  @property 
  def raw(self):
    """ Raw input for structure. """
    from ..error import ValueError
    from ..periodic_table import find as find_specie

    
    # symmgroup 
    result = str(self.symmgroup).rstrip().lstrip() + '\n'

    # number of atoms + atoms
    result += '{0}\n'.format(len(self.atoms))
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
    if not hasattr(value, '__iter__'): value = value.splitlines()
    self.symmgroup = int(value.pop(0))

    n = int(value.pop(0))
    self.atoms = []
    for line in value[:n]:
      line = line.split()
      type = int(line[0])
      if type < 100: type = find_specie(atomic_number=type).symbol
      self.add_atom(pos=[float(u) for u in line[1:4]], type=type)

  def read_input(self, tree, owner=None, **kwargs):
    """ Parses an input tree. """
    from ..tools.input import Tree
    from . import registered
    self[:] = []
    self.raw = tree.prefix
    kwargs['owner'] = self
    for key, value in tree:
      # parses sub-block.
      if isinstance(value, Tree): continue
      # check for special keywords.
      if key.lower() in registered:
        newobject = registered[key.lower()]()
      else: newobject = BaseKeyword(keyword=key)
      if len(value) > 0: 
        if hasattr(newobject, 'read_input'): 
          newobject.read_input(value, **kwargs)
        else:
          try: newobject.raw = getattr(value, 'raw', value)
          except: pass
      self.append(newobject)

  def append(self, keyword, raw=None):
    """ Appends an item.

        Some shortcuts are accpted:

        - 'breaksym'
        - 'keepsymm'
        - 'angstrom'
        - 'bohr'
        - 'fractional'
            
        :return: self, so calls can be chained. 
    """
    from . import registered
    if hasattr(keyword, 'lower')                                               \
       and keyword.lower() in [ 'breaksym', 'keepsymm', 'bohr', 
                                'fraction', 'angstrom' ]:
      keyword = registered[keyword]()
    super(Molecule, self).append(keyword, raw)
    return self
  def insert(self, i, keyword, raw=None):
    """ Inserts item at given position. 

        Some shortcuts are accpted:

        - 'breaksym'
        - 'keepsymm'
        - 'angstrom'
        - 'bohr'
        - 'fractional'
    """
    from . import registered
    if hasattr(keyword, 'lower')                                               \
       and keyword.lower() in [ 'breaksym', 'keepsymm', 'bohr', 
                                'fraction', 'angstrom' ]:
      keyword = registered[keyword]()
    super(Molecule, self).insert(i, keyword, raw)

  def __setitem__(self, i, value):
    """ Sets an item.

        Some shortcuts are accpted:

        - 'breaksym'
        - 'keepsymm'
        - 'angstrom'
        - 'bohr'
        - 'fractional'
    """
    from . import registered
    if hasattr(value, 'lower')                                                 \
       and value.lower() in [ 'breaksym', 'keepsymm', 'bohr', 
                              'fraction', 'angstrom' ]:
      value = registered[value](value=True)
    super(Molecule, self).__setitem__(i, value)

  @property
  def current_units(self):
    """ Current units in transformation list. 

        This will one of "bohr", "angstrom", "fractional", depending on the
        last encountered keyword in the list,
    """
    for u in self[::-1]:
      if u.keyword == 'fraction':   return 'fraction'
      elif u.keyword == 'bohr':     return 'bohr'
      elif u.keyword == 'angstrom': return 'angstrom'
    return 'angstrom'
  @property
  def is_fractional(self):
    """ True if current units are fractional. """
    return self.current_units == 'fraction'
  @property
  def is_bohr(self):
    """ True if current units are bohrs. """
    return self.current_units == 'bohr'
  @property
  def is_angstrom(self):
    """ True if current units are angstrom. """
    return self.current_units == 'angstrom'
  @property
  def is_breaksym(self):
    """ True if currently keeping symmetries. """
    for u in self[::-1]:
      if u.keyword == 'keepsymm':   return False
      elif u.keyword == 'breaksym': return True
    return True
    

  def eval(self):
    """ Evaluates current structure. 
    
        Runs crystal to evaluate current structure.

        :returns: a :py:class:`~pylada.crystal.cppwrappers.Structure` instance.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. error import GrepError
    from .. import crystal_program as program
    from .extract import Extract
    this = deepcopy(self)
    this.append('TESTGEOM')
    try:
      tmpdir = mkdtemp()

      with Changedir(tmpdir) as cwd:
        # writes input file
        with open('crystal.input', 'w') as file:
          file.write('{0}\n'.format(getattr(self, 'name', '')))
          file.write(this.print_input(filework=True))
        # creates process and run it.
        if hasattr(program, '__call__'): program = program()
        process = Process( program, stdin='crystal.input',
                           outdir=tmpdir, stdout='crystal.out',
                           stderr='crystal.err', dompi=False )

        process.start()
        process.wait()

        try: return Extract().input_structure
        except GrepError:
          message = this.print_input(filework=True) + '\n\n'
          with Extract().__stdout__() as file: message += file.read()
          raise ValueError(message + '\n\nInput structure is likely incorrect\n')
    finally:
      try: rmtree(tmpdir)
      except: pass

  @property
  def crystal_output(self):
    """ CRYSTAL program output. 
    
        This is a string which contains the CRYSTAL output from a geometry run.
        Electronic calculations are not performed.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. import crystal_program as program
    this = deepcopy(self)
    this.append('TESTGEOM')
    try:
      tmpdir = mkdtemp()

      with Changedir(tmpdir) as cwd:
        # writes input file
        with open('crystal.input', 'w') as file:
          file.write('{0}\n'.format(getattr(self, 'name', '')))
          file.write(this.print_input(filework=True))
        # creates process and run it.
        if hasattr(program, '__call__'): program = program(self)
        process = Process( program, stdin='crystal.input',
                           outdir=tmpdir, stdout='crystal.out',
                           stderr='crystal.err', dompi=False )

        process.start()
        process.wait()
        with open('crystal.out', 'r') as file: return file.read()
    finally:
      try: rmtree(tmpdir)
      except: pass

  @property
  def extprt(self):
    """ CRYSTAL external output.
    
        This is a string which contains the structure in CRYSTAL_'s external
        output format. Mostly for debugging.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. import crystal_program as program
    this = deepcopy(self)
    this.append('EXTPRT')
    this.append('TESTGEOM')
    try:
      tmpdir = mkdtemp()

      with Changedir(tmpdir) as cwd:
        # writes input file
        with open('crystal.input', 'w') as file:
          file.write('{0}\n'.format(getattr(self, 'name', '')))
          file.write(this.print_input(filework=True))
        # creates process and run it.
        if hasattr(program, '__call__'): program = program(self)
        process = Process( program, stdin='crystal.input',
                           outdir=tmpdir, stdout='crystal.out',
                           stderr='crystal.err', dompi=False )

        process.start()
        process.wait()
        with open('fort.34', 'r') as file: return file.read()
    finally:
      try: rmtree(tmpdir)
      except: pass

  @property
  def symmetry_operators(self):
    """ Symmetry operators, as determined by CRYSTAL.
    
        Runs crystal to evaluate the symmetry operators on the current
        structure.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. error import GrepError
    from .. import crystal_program as program
    from .extract import Extract
    this = deepcopy(self)
    this.append('TESTGEOM')
    try:
      tmpdir = mkdtemp()

      with Changedir(tmpdir) as cwd:
        # writes input file
        with open('crystal.input', 'w') as file:
          file.write('{0}\n'.format(getattr(self, 'name', '')))
          file.write(this.print_input(filework=True))
        # creates process and run it.
        if hasattr(program, '__call__'): program = program(self)
        process = Process( program, stdin='crystal.input',
                           outdir=tmpdir, stdout='crystal.out',
                           stderr='crystal.err', dompi=False )

        process.start()
        process.wait()

        try: return Extract().symmetry_operators
        except GrepError:
          message = this.print_input(filework=True) + '\n\n'
          with Extract().__stdout__() as file: message += file.read()
          raise ValueError(message + '\n\nInput structure is likely incorrect\n')
    finally:
      try: rmtree(tmpdir)
      except: pass

  def copy(self):
    """ Returns deep copy of self. """
    from copy import deepcopy
    return deepcopy(self)

  def print_input(self, **kwargs):
    """ Returns CRYSTAL input. """
    from warnings import warn
    from .input import print_input
    # corrects bug in crystal whereby mutliple atomsym commands will do
    # shit-loads of incorrect stuff, from overallocating array to operations
    # on unitialized arrays.
    map = self.output_map(**kwargs)
    indices = []
    for i, (key, value) in  enumerate(map[0][1]):
      if key.upper() == 'ATOMSYMM': indices.append(i)
    if len(indices) > 1:
      warn( "ATOMSYMM present more than once in input.\n"                       \
            "Multiple bugs in CRYSTAL make this unsafe.\n"                      \
            "All but the last instance of ATOMSYMM will be ignored.",
            UserWarning, 0 )
      indices.reverse()
      for i in indices[1:]: map[0][1].pop(i)
    # now print map.
    return print_input(map)
  def output_map(self, **kwargs):
    from ..tools.input import Tree
    from .input import find_sym, find_units
    if 'structure' not in kwargs: kwargs['structure'] = self
    if kwargs.get('breaksym', None) is None: kwargs['breaksym'] = True
    if kwargs.get('units', None) is None: kwargs['units'] = 'angstrom' 

    root = Tree()
    root[self.keyword] = Tree()
    result = root[self.keyword]
    if getattr(self, 'raw', None) is not None:
      result.prefix = str(self.raw)
    for item in self:
      keyword = getattr(item, 'keyword', None)
      value   = getattr(item, 'raw', None)
      if hasattr(item, 'output_map'):
        dummy = item.output_map(**kwargs)
        if dummy is not None:
          result.update(dummy)
          kwargs['breaksym'] = find_sym(result, result=kwargs['breaksym'])
          kwargs['units'] = find_units(result, result=kwargs['units'])
      elif value is not None:
        result[keyword] = value
      elif hasattr(value, '__iter__'):
        result[keyword] =  ' '.join(str(u) for u in value)
      else:
        result[keyword] = str(value)
      lowkey = keyword.lower()
      if lowkey == 'breaksym': kwargs['breaksym'] = True
      elif lowkey == 'keepsymm': kwargs['breaksym'] = False
      elif lowkey == 'angstrom': kwargs['units'] = 'angstrom'
      elif lowkey == 'bohr':     kwargs['units'] = 'bohr'
      elif lowkey == 'fraction': kwargs['units'] = 'fractional'
    return root
