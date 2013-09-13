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

""" Holds definition of GULP specific input classes. """
__docformat__ = "restructuredtext en"
__all__ = ['OptiKeyword', 'TwoBody', 'Species', 'Spring', 'Conv', 'Conp'] 
from collections import MutableMapping
from ..tools.input.keywords import BoolKeyword, BaseKeyword

class OptiKeyword(BoolKeyword):
  """ Keywords that depends upon opti for printing. """
  def output_map(self, **kwargs):
    """ Only prints if opti is on. """
    if kwargs['gulp'].opti:
      return super(OptiKeyword, self).output_map(**kwargs)
class Conv(OptiKeyword):
  """ Adds exclusive keyword behavior. """
  keyword = 'conv'
  """ Equivalent GULP keyword. """
  def __set__(self, instance, value):
    """ Makes sure only one of conv, conp is on. """
    if value is not None and value == True                                    \
       and getattr(instance, 'conp', False):
      instance.conp = False
    self.value = value
class Conp(OptiKeyword):
  """ Adds exclusive keyword behavior. """
  keyword = 'conp'
  """ Equivalent GULP keyword. """
  def __set__(self, instance, value):
    """ Makes sure only one of conv, conp is on. """
    if value is not None and value == True                                    \
       and getattr(instance, 'conv', False):
      instance.conv = False
    self.value = value
      

class Optimize(BoolKeyword):
  """ Sets optimization. """
  keyword = 'opti'
  """ Equivalent GULP keyword. """
  def output_map(self, **kwargs):
    """ Sets noflag keyword if no other opti keyword is set. """
    result = super(Optimize, self).output_map(**kwargs)
    if result is not None: 
      gulp = kwargs['gulp']
      if getattr(gulp, 'conv', None) is not True                               \
         and getattr(gulp, 'conp', None) is not True                           \
         and getattr(gulp, 'cellonly', None) is not True                       \
         and getattr(gulp, 'shell', None) is not True:
        result.update(noflags=True)  
    return result
                                        
class TwoBody(MutableMapping, BaseKeyword):
  """ Defines parameters for a two-body GULP functional. """
  def __init__(self, enabled=False, **kwargs):
    """ Creates the two-body functional. """
    super(TwoBody, self).__init__()
    self.params = {}
    """ Holds functional parameters. """
    self.enabled = enabled
    """ Whether this interaction is enabled or not. """
    
    # sets additional parameters to True or False
    for key, value in kwargs.iteritems(): setattr(self, key, value == True)

  def _regex_key(self, key):
    """ Returns regex-formatted value of a key. """
    from re import compile
    from ..error import IndexError
    regex = compile('([A-Z][a-z]?\d*)\s*\|?\s*(core|shell|)')
    groups = regex.match(key).groups()
    if groups is None:
      raise IndexError('Could not make sense of atom {0}.'.format(key))
    return groups[0] + '|' + ('core' if len(groups[1]) == 0 else groups[1])
  def _key(self, key):
    """ Creates a key from a two-body parameter. """
    if isinstance(key, str): return key
    if len(key) != 2: raise IndexError('Expected two atoms on input.')
    key = self._regex_key(key[0]), self._regex_key(key[1])
    return (key[1] + '-' + key[0]) if key[0] > key[1]                          \
           else (key[0] + '-' + key[1])
  def __getitem__(self, key):
    """ Retrieves two-body item. """
    return self.params[self._key(key)]
  def __setitem__(self, key, value):
    self.params[self._key(key)] = value
  def __delitem__(self, key, value):
    """ Deletes two-body interaction. """
    del self.params[self._key(key)] 
  def __len__(self):
    """ Number of independant two-body interactions. """
    return len(self.params)
  def __iter__(self):
    """ Iterates over the two-body interactions. """
    for k in self.params: yield k.split('-')

  def output_map(self, **kwargs):
    """ Returns output map. """
    # whether to print this or not.
    if self.enabled == False: return 
    # list of structures. 
    species = None
    if kwargs.get('structure', None) is not None:
      species = set([u.type for u in kwargs['structure']])
    # header is the keyword + additional parameters.
    header = self.keyword
    for key, value in self.__dict__.iteritems():
      if key[0] == '_' or key in ['params', 'keyword', 'enabled']: continue
      if not isinstance(value, bool): header += ' ' + str(value)
      elif value:                     header += ' ' + key
    # creates list of interactions.
    result = ""
    for key, value in self.params.iteritems():
      # if structure in input, check that interaction is needed.
      if species is not None:
        atom0, atom1 = key.split('-')
        atom0, atom1 = atom0.split('|')[0], atom1.split('|')[0]
        if atom0 not in species or atom1 not in species: continue
      # print interaction.
      result += key.replace('-', '   ').replace('|', ' ')  + '     '
      if hasattr(value, '__iter__'):
        result += ' '.join(str(float(u)) for u in value)
      else: result += str(float(value))
      result += '\n'
    return {header: result}

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly representation of the functional. """
    from ..tools.uirepr import template_ui_repr
    if exclude is None: exclude = ['params']
    else: exclude.append('params')

    results = template_ui_repr(self, imports, name, defaults, exclude)
    if defaults is None or self.enabled == True or len(self.params):
      results['{0}.enabled'.format(name)] = repr(self.enabled)
    for key, value in self.params.iteritems():
      index = ', '.join(repr(u) for u in key.split('-'))
      key = "{0}[{1}]".format(name, index)
      results[key] = str(float(value)) if not hasattr(value, '__iter__')       \
                     else ', '.join(str(u) for u in value) 
    return results

class Species(MutableMapping, BaseKeyword):
  """ Defines parameters for atomic parameters. """
  keyword = 'species'
  """ Equivalent GULP keyword. """
  def __init__(self, enabled=False, **kwargs):
    """ Creates the instance. """
    super(Species, self).__init__()
    self.params = {}
    """ Holds functional parameters. """
    self.enabled = enabled
    """ Whether this interaction is enabled or not. """
    
  def _key(self, key):
    """ Returns regex-formatted value of a key. """
    from re import compile
    from ..error import IndexError
    regex = compile('([A-Z][a-z]?\d*)\s*\|?\s*(core|shell|)')
    groups = regex.match(key).groups()
    if groups is None:
      raise IndexError('Could not make sense of atom {0}.'.format(key))
    return groups[0] + '|' + ('core' if len(groups[1]) == 0 else groups[1])
  def __getitem__(self, key):
    """ Retrieves atomic parameters for given specie. """
    return self.params[self._key(key)]
  def __setitem__(self, key, value):
    """ Sets atomic parameters for given specie. """
    self.params[self._key(key)] = value
  def __delitem__(self, key, value):
    """ Deletes atomic parameters for given specie. """
    del self.params[self._key(key)] 
  def __len__(self):
    """ Number of independant sets of atomic paramters. """
    return len(self.params)
  def __iter__(self):
    """ Iterates over atomic parameters. """
    for k in self.params: yield k.split('-')

  def output_map(self, **kwargs):
    """ Returns output map. """
    # whether to print this or not.
    if self.enabled == False: return 
    # list of structures. 
    species = None
    if kwargs.get('structure', None) is not None:
      species = set([u.type for u in kwargs['structure']])
    # header is the keyword + additional parameters.
    header = self.keyword
    for key, value in self.__dict__.iteritems():
      if key[0] == '_' or key in ['params', 'keyword', 'enabled']: continue
      if not isinstance(value, bool): header += ' ' + str(value)
      elif value:                     header += ' ' + key
    # creates list of interactions.
    result = ""
    for key, value in self.params.iteritems():
      # if structure in input, check that interaction is needed.
      if species is not None and key.split('|')[0] not in species: continue
      # print interaction.
      result += key.replace('|', ' ')  + '     '
      if hasattr(value, '__iter__'):
        result += ' '.join(str(float(u)) for u in value)
      else: result += str(float(value))
      result += '\n'
    return {header: result}

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly representation of the functional. """
    from ..tools.uirepr import template_ui_repr
    if exclude is None: exclude = ['params']
    else: exclude.append('params')

    results = template_ui_repr(self, imports, name, defaults, exclude)
    if defaults is None or self.enabled == True or len(self.params):
      results['{0}.enabled'.format(name)] = repr(self.enabled)
    for key, value in self.params.iteritems():
      index = ', '.join(repr(u) for u in key.split('-'))
      key = "{0}[{1}]".format(name, index)
      results[key] = str(float(value)) if not hasattr(value, '__iter__')       \
                     else ', '.join(str(u) for u in value) 
    return results

class Springs(Species):
  keyword = 'spring'
  """ GULP equivalent keyword. """
  def _key(self, key):
    """ Returns regex-formatted value of a key. """
    from re import compile
    from ..error import IndexError
    regex = compile('([A-Z][a-z]?\d*)')
    groups = regex.match(key).groups()
    if groups is None:
      raise IndexError('Could not make sense of atom {0}.'.format(key))
    return groups[0] 
