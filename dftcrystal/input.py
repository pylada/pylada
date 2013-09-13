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
__all__ = ['AttrBlock']
from ..tools.input import AttrBlock as AttrBlockBase, BaseKeyword,             \
                          BoolKeyword as BaseBoolKeyword

class AttrBlock(AttrBlockBase):
  """ Defines block input to CRYSTAL. 
  
      A block is a any set of input which starts with a keyword and ends with
      an END. It can contain other sub-keywords.
      This particular flavor creates attributes from the  inner input keywords.
      This supposes that each keyword is only ever inputed once. 

      It can contain subitems.
  """
  def __init__(self, keyword=None, raw=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    super(AttrBlock, self).__init__(keyword=keyword, raw=raw)
    
  def output_map(self, **kwargs):
    """ CRYSTAL input.
    
        The input is composed as followed:
        
          - the keyword of this block, if it exists.
          - the 'raw' input, if is exists.
          - each input keyword to CRYSTAL contained by the block.
    """
    from ..tools.input import Tree
    raw = getattr(self, 'raw', None)
    result = super(AttrBlock, self).output_map(**kwargs)
    if result is None:
      if raw is None: return None
      result = Tree()
    if raw is not None: result.prefix = raw
    return result

  def read_input(self, tree, owner=None, **kwargs):
    """ Parses an input tree. """
    from ..error import internal
    from ..tools.input import Tree
    from .input import find_sym, find_units
    from . import registered
    if hasattr(self, 'raw') and len(getattr(tree, 'prefix', '')) > 0: 
      try: self.raw = tree.prefix
      except: pass
    breaksym = kwargs.get('breaksym', False)
    units    = kwargs.get('units', 'angstrom')
    for key, value in tree:
      # parses sub-block.
      # check for special keywords.
      lowkey = key.lower()
      if lowkey in self._input: 
        newobject = self._input[lowkey] 
        if hasattr(newobject, 'read_input'):
          newobject.read_input( value, owner=self,
                                breaksym=breaksym,
                                units=units )
        elif hasattr(newobject, 'raw'): newobject.raw = value
        elif hasattr(newobject, '__set__'): newobject.__set__(self, value)
        else:
          raise internal( "Pylada doesn't understand how to read input to {0}"   \
                          .format(lowkey) )
        continue
      if lowkey == 'breaksym': breaksym = True; continue
      if lowkey == 'keepsymm': breaksym = False; continue
      if lowkey == 'bohr':       units = 'bohr'; continue
      if lowkey == 'angstrom':   units = 'angstrom'; continue
      if lowkey == 'fractional': units = 'fractional'; continue
      if isinstance(value, Tree):
        newobject = registered.get(lowkey, AttrBlock)()
        newobject.read_input(value, owner=self, breaksym=breaksym, units=units)
        breaksym = find_sym(value, result=breaksym)
        units = find_units(value, result=units)
      # creates new object.
      elif lowkey in registered:
        newobject = registered[lowkey]()
        if len(value) > 0:
          try: newobject.raw = getattr(value, 'raw', value)
          except: pass
      elif len(value) == 0: newobject = True
      else: newobject = value
      self.add_keyword(lowkey, newobject)

  def _read_nested_group(self, tree, owner=None, **kwargs):
    """ Creates nested groups.

        It is not easy to guess the type of a yet to be created object.
        This method is used by derived classes to create default nested groups.
    """
    result = self.__class__()
    result.read_input(tree, owner=self, **kwargs)
    return result


class BoolKeyword(BaseBoolKeyword):
  """ Bool option for CRYSTAL.
  
      Makes sure it is set to True from its mere presence in the input file.
  """
  def read_input(self, *args, **kwargs): self.value = True


class SetPrint(BaseKeyword):
  """ Print options. """
  keyword = "setprint"
  def __init__(self, args=None):
    super(SetPrint, self).__init__()
    self.options = args.copy() if args is not None else {}
  def __getitem__(self, name):
    return self.options[name]
  def __setitem__(self, name, value):
    self.options[name] = value
  def output_map(self, **kwargs):
     if len(self.options) == 0: return None
     result = str(len(self.options)) + '\n'
     for i, (key, value) in enumerate(self.options.iteritems()):
       if i != 0 and  i % 5 == 0: result += '\n'
       result += "{0} {1}  ".format(key, value)
     return {self.keyword: result}
  def __get__(self, instance, owner=None): return self
  def __set__(self, instance, value):
    self.options = {} if value is None else value.copy()
  def __delete__(self, instance): self.options = {}
  def __len__(self): return len(self.options)

  def read_input(self, value, **kwargs):
    from ..error import ValueError
    self.options = {}
    if value is None: return
    value = value.split()
    if len(value) < 3: return
    N = int(value[0]) * 2 + 1
    if len(value) < N:
      raise ValueError('Array of options is too small in SetPrint.')
    for key, value in zip(value[1:N:2], value[2:N:2]):
      self.options[int(key)] = int(value)
  def __repr__(self):
    """ Dumps representation to string. """
    if len(self.options) == 0: return "{0.__class__.__name__}()".format(self)
    return "{0.__class__.__name__}({0.options!r})".format(self)

def print_input(map):
  """ Prints output map to string. """
  from ..tools.input import Tree
  result = ''
  for key, item in map.iteritems():
    if isinstance(item, Tree):
      prefix = getattr(item, 'prefix', None)
      if hasattr(prefix, '__len__') and len(prefix) == 0: prefix = None
      if key != 'BASISSET': result += key.upper()  + '\n'
      if prefix is not None: result += prefix.upper().rstrip() + '\n'
      dummy = print_input(item).rstrip().lstrip()
      if len(dummy) > 0: result += dummy.upper() + '\n'
      result += 'END ' + key.upper() + '\n'
    else: 
      if item in ['True', 'None']: 
        result += str(key).upper().rstrip().lstrip() + '\n'
      elif item is False: continue
      elif item is True or item is None:
        result += key.upper() + '\n'
      elif item is not None and len(item) > 0 and item != 'False':
        result += key.upper().rstrip().lstrip() + '\n'                         \
                  + item.upper().rstrip() + '\n'
  return result

def find_units(tree, breakpoint=None, result=None):
  """ Figures out current units in tree stream. """
  from ..tools.input import Tree
  if breakpoint is not None: breakpoint = breakpoint.lower()
  for key, value in tree:
    key = key.lower()
    if breakpoint is not None and key == breakpoint: return result
    elif key == 'fraction': result = 'fraction'
    elif key == 'bohr':       result = 'bohr'
    elif key == 'angstrom':   result = 'angstrom'
    elif isinstance(value, Tree): result = find_units(value, result=result)
  return result

def find_sym(tree, breakpoint=None, result=None):
  """ Figures out current units in tree stream. """
  from ..tools.input import Tree
  if breakpoint is not None: breakpoint = breakpoint.lower()
  for key, value in tree:
    key = key.lower()
    if breakpoint is not None and key == breakpoint: return result
    elif key == 'keepsymm': result = False
    elif key == 'breaksym': result = True
    elif isinstance(value, Tree): result = find_sym(value, result=result)
  return result
