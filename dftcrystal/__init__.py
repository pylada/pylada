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

""" Python wrapper for the CRYSTAL dft code. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Extract', 'AffineTransform', 'DisplaceAtoms', 'InsertAtoms',
            'ModifySymmetry', 'RemoveAtoms', 'Slabcut', 'read',
            'Slabinfo','Crystal', 'Molecule', 'Slab', 'Functional', 'Shell',
            'read_gaussian_basisset', 'MassExtract', 'read_input', 'exec_input',
            'Elastic', 'External', 'Supercell', 'Supercon', 'KeepSymm', 'BreakSym',
            'BohrUnits', 'AngstromUnits', 'Fractional' ]

from .basis import Shell
from .functional import Functional
from .extract import Extract, MassExtract
from .geometry import AffineTransform, DisplaceAtoms, InsertAtoms,              \
                      ModifySymmetry, RemoveAtoms, Slabcut, Slabinfo, Slab,     \
                      Elastic, Supercell, Supercon, BreakSym, KeepSymm,         \
                      BohrUnits, AngstromUnits, Fractional
from .crystal import Crystal
from .external import External
from .molecule import Molecule
registered = { 'atomrot':    AffineTransform,
               'atomdisp':   DisplaceAtoms,
               'atominse':   InsertAtoms,
               'modisymm':   ModifySymmetry,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'elastic':    Elastic,
               'supercel':   Supercell,
               'supercon':   Supercon,
               'crystal':    Crystal,
               'external':   External,
               'molecule':   Molecule,
               'breaksym':   BreakSym,
               'keepsymm':   KeepSymm,
               'bohr':       BohrUnits,
               'angstrom':   AngstromUnits,
               'fraction':   Fractional }
""" Keywords used in creating the geometry. """

def read(path):
  """ Reads CRYSTAL input from file

      Reads a file and finds the CRYSTAL_ input within. It then parses the file
      and creates a :py:class:`~functional.Functional` and
      :py:class:`~molecule.Molecule` (or derived) object.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      :returns: A two-tuple containing the structure and the functional

      :rtype: (:py:class:`~molecule.Molecule`, :py:class:`~functional.Functional`)

      .. _CRYSTAL: http://www.crystal.unito.it/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from .parse import parse
  b = parse(path)
  key = b.keys()[0]

  func = Functional()
  crys = Crystal()
  func.read_input(b)
  crys.read_input(b[key]['CRYSTAL'])
  crys.name = key
  if func.optgeom.enabled and len(crys) > 0                                    \
     and crys[-1].keyword.lower() in ['breaksym', 'keepsymm']:
     crys.pop(-1)
  return crys, func


def read_gaussian_basisset(path):
  """ Reads a GAUSSIAN94 basis set

      This is meant to read input files from the `basis set exchange website`__
      in the GAUSSIAN94_ format.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      .. __: https://bse.pnl.gov/bse/portal
      .. _GAUSSIAN94: http://www.gaussian.com/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from ..error import GrepError
  from ..misc import RelativePath
  if isinstance(path, str): 
    if path.find('\n') == -1:
      with open(RelativePath(path).path) as file: return read_gaussian_basisset(file)
    else:
      return read_gaussian_basisset(path.split('\n').__iter__())

  for line in path: 
    if set(['*']) == set(line.rstrip().lstrip()): break

  # read specie type
  try: line = path.next().split()
  except StopIteration: raise GrepError('Unexpected end of file')
  specie = line[0]
  result = {specie: []}
  
  try: 
    while True:
      line = path.next().split()
      if len(line) != 3:
        line = path.next().split()
        if len(line) != 2: break
        specie = line[0]
        result[specie] = []
        continue
      type, n, scale = line[0].lower(), int(line[1]), float(line[2])
      shell = Shell(type)
      for i in xrange(n): shell.append(*path.next().split())
      result[specie].append(shell)

  except StopIteration: pass

  return result

def read_input(filepath="input.py", namespace=None):
  """ Specialized read_input function for dftcrystal. 
  
      :Parameters: 
        filepath : str
          A path to the input file.
        namespace : dict
          Additional names to include in the local namespace when evaluating
          the input file.

      It add a few names to the input-file's namespace. 
  """
  from ..misc import read_input

  # names we need to create input.
  input_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input':
      input_dict[k] = globals()[k]
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)

def exec_input( script, global_dict=None, local_dict=None,
                paths=None, name=None ):
  """ Specialized exec_input function for vasp. """
  from ..misc import exec_input

  # names we need to create input.
  if global_dict is None: global_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input':
      global_dict[k] = globals()[k]
  return exec_input(script, global_dict, local_dict, paths, name)

from pylada import is_interactive
if is_interactive:
  from .interactive import *

del is_interactive
