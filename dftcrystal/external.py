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
__all__ = ['External']
from .molecule import Molecule
class External(Molecule):
  """ Functional CRYSTAL-wise structure, starting from Pylada structure.

      Provides a mixture of CRYSTAL_ and Pylada-style structures. The starting
      structure is a :py:class:`~pylada.crystal.cppwrappers.Structure` instance
      which is inputed into CRYSTAL_ via the EXTERNAL keyword. This instance
      also accepts functional modifications which will act upon the initial
      structure, much as :py:class:`~pylada.dftcrystal.crystal.Crystal`
      instances.

      There are two ways of initializing an instance. The first allows the same
      interface as :py:class:`~pylada.crystal.cppwrappers.Structure`.

      >>> from quantities import angstrom
      >>> from pylada.dftcrystal import External
      >>> external = External([[0, 0.5, 0.5],
      ...                      [0.5, 0, 0.5],
      ...                      [0.5, 0.5, 0], scale=5.45*angstrom)            \\
      ...                    .add_atom(0,0,0, 'Si')                           \\
      ...                    .add_atom(0.25, 0.25, 0.25, 'Si')

      The second approach is to use a predefined
      :py:class:`~pylada.crystal.cppwrappers.Structure` instance:

      >>> from pylada.dftcrystal import External
      >>> external = External(copy=diamond)

      Where ``diamond`` in the snippet above is the instance in question. Note
      that it is deepcopied upon initialization.  In both cases, the initial
      structure can be accessed as :py:attr:`initial`. 

      In both cases, CRYSTAL_ transformations can be added to modify the
      initial structure.

      >>> external.append('bohr')
      >>> external.append(DiplaceAtoms().add_atom(0.1, 0, 0, 1))

      In the snippet above, the displacement is in the current units of the
      structure, as understood by CRYSTAL_. E.g. in "bohr"'s here. In other
      words, the transformations should reference the structure the way
      CRYSTAL_ understands it. It might be indicated to use :py:meth:`eval`.

      .. note::
         
         Pylada will attempt to discover the symmetries at run time. It is also
         possible to enter them explicitely by setting :py:attr:`initial`'s
         ``spacegroup`` attribute to a list of 4x3 matrices.
  """
  keyword = 'external'
  """ CRYSTAL keyword. """
  def __init__(self, *args, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    from ..crystal import Structure
    self.__dict__['initial'] = None
    super(External, self).__init__()
    del self.symmgroup
    del self.__dict__['atoms']

    self.initial = kwargs.pop('copy', None)
    """ Initial structure. 
    
        This is a :py:class:`~pylada.crystal.cppwrappers.Structure` instance.
    """
    if self.initial is None: self.initial = Structure(*args, **kwargs)
    else: self.initial = self.initial.copy()

  def add_atoms(self, *args, **kwargs): 
    """ Adds atom to initial structure. """
    self.initial.add_atom(*args, **kwargs)
    return self
  @property
  def atoms(self):
    """ Alias to the initial structure. """
    return self.initial
  @atoms.setter
  def atoms(self, value):
    """ Alias to the initial structure. """
    from ..error import AttributeError
    if self.initial is not None: 
      raise AttributeError('Cannot set atoms in External instance.')
    self.__dict__['atoms'] = value

  @property
  def cell(self):
    """ Cell-vectors of the initial instance. """
    return self.initial.cell
  @cell.setter
  def cell(self, value):
    self.initial.cell = value

  def __repr_header__(self):
    """ Dumps representation to string. """
    # prints construction part.
    args = ""
    if self.initial is not None:
      copy = self.initial.copy()
      copy.clear()
      args = repr(copy)
      args = args[args.find('(')+1:args.rfind(')')]
    return '{0.__class__.__name__}({1})'.format(self, args)

  @property 
  def raw(self):
    """ Raw input for structure. """
    return ""
  @raw.setter 
  def raw(self, value): pass

  def output_map(self, **kwargs):
    """ Writes external file and calls base class. """
    from ..crystal import write
    if kwargs.get('filework', False) == True:  write.crystal(self.initial)
    return super(External, self).output_map(**kwargs)
