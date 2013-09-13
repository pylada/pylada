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

""" Methods to parse CRYSTAL input. """
__docformat__ = "restructuredtext en"
__all__ = ['parse']

def parse(path, needstarters=True):
  """ Reads crystal input. """
  from re import compile
  from ..misc import RelativePath
  from ..error import IOError
  from ..tools.input import Tree
  from .. import CRYSTAL_input_blocks as blocks, CRYSTAL_geom_blocks as starters
  if isinstance(path, str): 
    if path.find('\n') == -1:
      with open(RelativePath(path).path) as file: return parse(file,needstarters)
    else:
      return parse(path.split('\n').__iter__(),needstarters)

  exc = [ 'LYP', 'P86', 'PBE', 'PBESOL', 'PWLSD', 'PZ', 'VBH', 'WL', 'VWN',
          'BECKE', 'LDA', 'PWGGA', 'SOGGA', 'WCGGA', 'INPUT' ]
  keyword_re = compile('^(?:[A-Z](?!\s))')

  if needstarters: 
    title = None
    for i, line in enumerate(path):
      if len(line.split()) > 0 and line.split()[0] in starters: break
      title = line
    if title is None: raise IOError('Could not find CRYSTAL input.')
    if len(line) == 0: raise IOError('Could not find CRYSTAL input.')
    elif line.split()[0] not in starters: 
      raise IOError('Could not find CRYSTAL input.')
    title = title.rstrip().lstrip()
    nesting = [title, line.split()[0]]
    results = Tree()
    keyword = line.split()[0]
    raw = ''
  else: 
    nesting = []
    results = Tree()
    keyword = None
    raw = ''
  # reads crystal input.
  for line in path:
    # special case of INPUT keyword
    if len(line.split()) == 0: continue
    # found keyword
    elif keyword_re.match(line) is not None and line.split()[0] not in exc:
      newkeyword = line.split()[0]
      current = results.descend(*nesting)

      # first of subblock
      if len(nesting) and keyword == nesting[-1]: current.prefix = raw
      elif keyword is not None and keyword[:3] != 'END'                        \
           and keyword[:6] != 'STOP':
        current.append((keyword, raw)) 

      # normal keyword
      if newkeyword in blocks                                                  \
         and not (newkeyword == 'SLAB' and nesting[-1] == 'CRYSTAL'): 
        nesting.append(newkeyword)
      # found end, pop nesting.
      if newkeyword[:3] == 'END' or newkeyword[:6] == 'STOP': 
        current = nesting.pop(-1)
        if current in starters:
          nesting.append('BASISSET')
          newkeyword = 'BASISSET'
        elif current == 'BASISSET': 
          nesting.append('HAMSCF')
          newkeyword = 'HAMSCF'
        elif current == 'HAMSCF': current = nesting.pop(-1)
        if len(nesting) == 0: break
      raw = ''
      keyword = newkeyword
    # found raw string
    else:
      raw += line.rstrip().lstrip()
      if raw[-1] != '\n': raw += '\n'
  return results
