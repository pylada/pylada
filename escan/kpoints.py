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

""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ['KPoints', 'KGrid', 'ReducedKGrid', 'ReducedKDensity']
from abc import ABCMeta, abstractmethod
 
class KPoints(object):
  """ Abstract base class for callable KMesh objects. 

      It is expected that the kpoints are expressed in cartesian coordinates
      of the reciprocal space. There is no 2|pi|. In other words, the
      following code holds true:

      >>> for count, kvec in kpoints_instance(structure): 
      >>>   assert abs(exp(-2e0*pi*1j*dot(kvec, structure.cell[:,0])) -1e0) < 1e-12

      The code above assumes that structure is a valid `crystal.Structure`
      object.

      .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
  """
  __metaclass__ = ABCMeta

  def __init__(self): super(KPoints, self).__init__()

  @abstractmethod
  def _mnk(self, input, output):
    """ Returns iterable yielding (multiplicity, kpoints). 
    
        :Parameters:
          input : `pylada.crystal.Structure`
            The structure before vff relaxation (if any). 
          output : `pylada.crystal.Structure`
            The structure after vff relaxation. If no relaxation was performed,
            then should be the same structure as the input structure.
    """
    pass
  @abstractmethod
  def __repr__(self):
    """ This object must be representable. """
    pass
  
  def multiplicity(self, input, output):
    """ Generator yielding multiplicity. """
    for count, k in self._mnk(input, output): yield count
  def kpoint(self, input, output):
    """ Generator yielding kpoints. """
    for count, k in self._mnk(input, output): yield k
  def __call__(self, input, output): 
    """ Iterator over (multiplicity, kpoint) tuples. """
    for r in self._mnk(input, output): yield r
    

class KContainer(object):
  """ Simple KPoints class which acts as a container. """
  def __init__(self, kpoints, multiplicity, relax=True):
    """ Initializes the kpoint container. """
    self.kpoints = [k for k in kpoints]
    """ Sequence of kpoints. """
    self.multiplicity = multiplicity
    """ Sequence with the multiplicity of the respective kpoints. """
    if self.multiplicity is None:
      self.multiplicity = [1e0 / len(self.kpoints) for k in self.kpoints]
    else: self.multiplicity = [m for m in self.multiplicity]
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  def _mnk(self, input, output):
    """ Loop over array of kpoints. 
    
        If ``self.relax`` is True, then kpoints are relaxed from the input cell
        to the relaxed output cell.
    """
    from numpy import dot
    from numpy.linalg import inv

    if self.relax: deformation = dot(inv(output.cell.T), input.cell.T)
    for count, k in zip(self.kpoints, self.multiplicity):
      if self.relax: k = dot(deformation, k)
      yield count, k
    
  def kpoint(self, input, output):
    """ Generator yielding kpoints. """
    for count, k in self._mnk(input, output): yield k
  def __call__(self, input, output): 
    """ Iterator over (multiplicity, kpoint) tuples. """
    for r in self._mnk(input, output): yield r

  def __repr__(self, *args):
    return '{0.__class__.__name__}({1},{2})'\
           .format(self, repr(self.kpoints), repr(self.multiplicity))

class KGrid(KPoints):
  """ Unreduces kpoint grid with offsets. """

  def __init__(self, grid = None, offset = None, relax = True):
    """ Initializes unreduced KGrid. """
    from numpy import array
    KPoints.__init__(self)
    self.grid = grid if grid is not None else array([1,1,1])
    """ Grid dimensions in reciprocal space. """
    self.offset = offset if offset is not None else array([0,0,0])
    """ Offset from Gamma of the grid. """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """


  def _mnk(self, input, output):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv
    from numpy import zeros, array, dot, multiply
    from ..crystal.gruber import Reduction
    
    invgrid = array([1./float(u) for u in self.grid])
    offset = multiply(array(self.offset), invgrid) - array([0.5,0.5,0.5])
    cell = inv(Reduction()(output.cell, recip=True)).T
    if self.relax:
      deformation = dot(output.cell, input.cell.T)
      cell = dot(deformation, cell)

    a = zeros((3,), dtype='float64')
    weight = 1e0 / float(self.grid[0] * self.grid[1] * self.grid[2])
    for x in xrange(self.grid[0]):
      for y in xrange(self.grid[1]):
        for z in xrange(self.grid[2]):
          a = dot(cell, multiply(array([x,y,z]), invgrid) + offset)
          yield weight, array(a.flat)
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_one = all( abs(array(self.grid)-array([1,1,1])) < 1e-12 )
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_one and is_zero:
      return '{0.__class__.__name__}(relax={0.relax})'.format(self)
    if is_one:
      return '{0.__class__.__name__}(relax={0.relax}, '\
             'offset=({0.offset[0]},{0.offset[0]},{0.offset[0]}))'\
             .format(self)
    if is_zero:
      return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}), '\
             'relax={0.relax}))'.format(self)
    return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}), '\
           '({0.offset[0]},{0.offset[0]},{0.offset[0]}), relax={0.relax})'.format(self)


class KDensity(KGrid):
  """ Unreduced kpoint grid parameterized by the density and offset. """

  def __init__(self, density, offset = None, relax=True):
    """ Initializes unreduced KGrid. """
    KGrid.__init__(self, relax=relax, offset=offset)
    self.density = density
    """ 1-dimensional density in cartesian coordinates (1/Angstrom). """

  def _mnk(self, input, output):
    """ Yields kpoints on the grid. """
    from numpy import floor
    from numpy.linalg import inv, norm, det
    from ..crystal.gruber import Reduction
    
    cell = inv(Reduction()(output.cell,recip=True) * output.scale).T
    assert abs(det(inv(output.cell*output.scale)) - det(cell)) < 1e-8
    density = self.density
    if hasattr(density, 'rescale'): density.rescale(1e0/Angstrom)

    self.grid = [int(max(1, floor(norm(a) * self.density+0.5))) for a in cell]
    result = KGrid._mnk(self, input, output)
    return result
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_zero:
      return '{0.__class__.__name__}({0.density}, relax={0.relax})'.format(self, repr(self.cell))
    return '{0.__class__.__name__}({1}, ({2[0]},{2[0]},{2[0]}), relax={0.relax})'\
           .format(self, self.density, self.offset)

def _reduced_grids_factory(name, base):
  class ReducedKGrid(base): 
    def __init__(self, *args, **kwargs):
      """ Initializes reduces k-grid.
      
          :param args: Passed on to base class.
          :param kwargs: Passed on to base class.
          :kwarg tolerance: real number which defines the criterion by which
            two vectors are recognized as equal. Defaults to 1e-12.
      """
      self.tolerance = kwargs.pop('tolerance', 1e-12)
      """ Criteria to determine whether two k-vectors are the same. """
      base.__init__(self, *args, **kwargs)

    def iter_equivalents(self, input, output):
      """ Yields iterators over equivalent kpoints.
      
          :Parameters:
            input 
              Input structure as originally given before structural relaxation.
            output
              Output structure after structural relaxation, as used in
              electronic structure calculation.

          :return:
            Yields an iterator over equivalent kpoints 
            The iterator themselves yield 4-tuples:

            - index in unreduced mesh
            - multiplicity of unreduced kpoints
            - unreduced kpoint
            - operation to go to reduced kpoint. 
      """
      from numpy import dot
      from numpy.linalg import inv
      from ..crystal import to_origin, to_voronoi, SymmetryOperator
      lattice = output.to_lattice()
      kcell = inv(lattice.cell).T

      # now checks whether symmetry kpoint exists or not.
      seen = []
      for j, (mult, kpoint) in enumerate(base._mnk(self, input, output)): 
        found = False
        kpoint = to_voronoi(kpoint, kcell)
        for i, others in enumerate(seen):
          index, m, vec, op = others[0]
          for op in lattice.space_group:
            try: 
              u = to_origin(dot(op.op, kpoint), kcell, vec)
            except ValueError:
              print kpoint
              raise
            if all(abs(u) < self.tolerance):
              found = True
              seen[i].append((j, mult, kpoint, op))
              break
          if found: break
        if found == False: seen.append([(j, mult, kpoint.copy(), SymmetryOperator())])

      for i, kpoints in enumerate(seen): yield kpoints.__iter__()

    def _mnk(self, input, output):
      """ Returns list of inequivalent vectors with multiplicity. """
      for equivs in self.iter_equivalents(input, output):
        index, m0, k0, op = equivs.next()
        for index, m, k, op in equivs: m0 += m
        yield m0, k0

    def unreduced(self, input, output):
      """ Yields unreduced kpoints. """
      for mult, kpoint in base._mnk(self, input, output): yield mult, kpoint

    def mapping(self, input, output):
      """ Yields index of unreduced kpoint in array of reduced kpoints. """
      from operator import itemgetter
      indices = []
      for i, equivs in enumerate(self.iter_equivalents(input, output)):
        indices.extend([(i, index) for index, m, k, op in equivs])
      for i, index in sorted(indices, key=itemgetter(1)): yield i

    def __repr__(self):
      """ Represents this object. """
      if self.tolerance == 1e-12: return base.__repr__(self)
      result = base.__repr__(self)
      result = result[:-1].rstrip()
      if result[-1] == '(': return result + 'tolerance={0})'.format(self.tolerance)
      return result + ', tolerance={0})'.format(self.tolerance)

  ReducedKGrid.__name__ = name
  ReducedKGrid.__doc__ = """ {0} reduced according to symmetries. """.format(base.__name__)
  ReducedKGrid.__module__ = base.__module__ 
  return ReducedKGrid


ReducedKGrid    = _reduced_grids_factory('ReducedKGrid', KGrid)
ReducedKDensity = _reduced_grids_factory('ReducedKDensity', KDensity)
