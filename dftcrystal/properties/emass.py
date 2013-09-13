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

""" Computes effective masses for any structure. """
from ...tools.makeclass import makeclass, makefunc
from ...tools import make_cached
from .extract import Extract as ExtractBase
from . import Properties

class Extract(ExtractBase):
  """ Extraction for effective mass. """
  def __init__(self, *args, **kwargs):
     super(Extract, self).__init__(*args, **kwargs)

  @property 
  def success(self):
    """ True if successful run. 
    
        Checks this is an effective mass calculation.
    """
    try: self._details
    except: return False
    return ExtractBase.success.__get__(self)

  @property
  @make_cached
  def _details(self):
    """ Parameters when calling the effective mass routine. """
    from re import compile
    from ...misc import exec_input
    from ...error import GrepError
    start = compile(r'^#+ EMASS DETAILS #+$')
    end = compile(r'^#+ END EMASS DETAILS #')
    with self.__stdout__() as file:
      lines = None
      for line in file:
        if start.match(line): lines = ""; break
      if lines is None: raise GrepError('Could not find call parameters.')
      for line in file:
        if end.match(line):  break
        lines += line

    input = exec_input(lines)
    return { 'center': input.center,
             'nbpoints': input.nbpoints,
             'range': input.range }
  @property
  def center(self): return self._details['center']
  @property
  def nbpoints(self): return self._details['nbpoints']
  @property
  def range(self): return self._details['range']
  @property
  def kpoints(self): return self.bandstructure.kpoints
  @property
  def eigenvalues(self): return self.bandstructure.eigenvalues

  @staticmethod
  def _orders(orders):
    """ Order up to which taylor coefficients should be computed. """
    result = orders 
    if result is None: result = [0, 2]
    if not hasattr(result, '__iter__'): result = [result]
    return sorted(result)

  def parameters(self, orders):
    """ Parameters for the least-square-fit approach. """
    from numpy import zeros, pi, dot
    from numpy.linalg import inv
    from math import factorial

    # standardize order input. Now a list.
    order = self._orders(orders)
    # computes number of free parameters, accounting for analyticity.
    u = 0
    for order in orders:
      for p in xrange(order+1):
        for j in xrange(p+1): u += 1
    # creates parameter array
    result = zeros((len(self.kpoints), u), dtype="float64")
    # populates it.
    u = 0
    recipcell = inv(self.structure.cell).T
    kpoints = dot(recipcell, (self.kpoints - self.center).T).T * 2.0 * pi      \
              / self.structure.scale 
    stepsize = zeros(len(self.kpoints), dtype="float64")
    for order in orders:
      for p in xrange(order+1):
        for q in xrange(p+1):
          # constructs array of (x_j-x_0)^q_i, with multiplication over i
          # representing the cartesian coordinates and j the kpoints.
          stepsize[:]                                                          \
            = 1. / (factorial(order-p)*factorial(p-q)*factorial(q))
          for i, j in zip([order-p, p-q, q], xrange(3)): 
            if i != 0: stepsize *= kpoints[:,j]**i
          # fill in the parameters for this expansion.
          result[:, u] = stepsize 
          u += 1
    return result

  def fit(self, orders=None):
    """ Returns all fitted taylor coefficients. """
    from numpy import sum, argmin
    from numpy.linalg import lstsq
    if orders is None: orders = [0, 2]
    # computes measurements
    eigenvalues = self.eigenvalues
    if 0 not in orders: # recenter.
      index = argmin(sum((self.kpoints - self.center)**2, axis=1))
      if len(eigenvalues.shape) == 3: # spin  polarized
        eigenvalues[0] -= eigenvalues[0, index]
      else: # spin unpolarized
        eigenvalues -= eigenvalues[index]
    return lstsq(self.parameters(orders), eigenvalues.magnitude)

  def tensors(self, orders=None):
    """ Return tensors for each band. """
    from numpy import zeros, array
    from itertools import chain, permutations
    from quantities import angstrom
    orders = self._orders(orders)
    all_xs = self.fit(orders)[0].T
    result = [[]]*(len(orders))
    
    current_index, current_order = 0, 0
    # add results for order 0.
    if 0 in orders:
      result[current_order] = all_xs[:, current_index].copy()                  \
                              * self.eigenvalues.units                         
      current_index += 1                                                       
      current_order += 1                                                       
    # add results for order 1.                                                 
    if 1 in orders:                                                        
      result[current_order] = all_xs[:, current_index:current_index+3]         \
                              * self.eigenvalues.units * angstrom
      current_index += 3
      current_order += 1

    # compute index ranges for subsequent orders.
    u, indices_range = current_index, [current_index]
    for order in orders[current_order:]:
      for p in xrange(order+1):
        for q in xrange(p+1): u += 1
      indices_range.append(u)

    # loop over remaining orders.
    for order, startu in zip(orders[current_order:], indices_range):
      # loop over results for each band.
      for band_x in all_xs:
        u = startu
        # create tensor from n vector.
        dummy = zeros([3]*order, dtype="float64")
        for p in xrange(order+1):
          for q in xrange(p+1):
            indices = [ [i]*j for i, j in
                         zip([0, 1, 2], [order-p, p-q, q]) if j != 0 ]
            indices = set(permutations(chain(*indices)))
            for index in indices:
              dummy[index] = band_x[u] if abs(band_x[u]) > 1e-12 else 0
            u += 1
        # add tensor to results for that order.
        result[current_order].append(dummy)
      result[current_order] = array(result[current_order])                     \
                              * self.eigenvalues.units * angstrom**order
      current_order += 1

    # got it all.
    return result

  @property
  def breakpoints(self):
    """ Indices for start of each path. """
    from numpy import any, abs, cross
    breakpoints, last_dir = [0], None
    for i, k in enumerate(self.kpoints[1:]):
      if last_dir is None: last_dir = k - self.kpoints[breakpoints[-1]]
      elif any( abs(cross(last_dir, k-self.kpoints[breakpoints[-1]])) > 1e-8):
        breakpoints.append(i+1)
        last_dir = None
    return breakpoints + [len(self.kpoints)]

  @property
  def directions(self):
    """ Direction for each path. """
    from numpy import dot, array
    from numpy.linalg import norm, inv
    from quantities import angstrom
    results = []
    breakpoints = self.breakpoints
    recipcell = inv(self.structure.cell).T
    for start, end in zip(breakpoints[:-1], breakpoints[1:]):
      results.append(self.kpoints[end-1] - self.kpoints[start])
      results[-1] = dot(recipcell, results[-1])
      results[-1] /= norm(results[-1])
    return array(results) / angstrom

  def emass(self, orders=None):
    """ Computes effective mass for each direction. """
    from numpy import dot, concatenate, pi, array
    from numpy.linalg import inv, lstsq
    from math import factorial
    from quantities import angstrom, emass, h_bar
    from ...error import ValueError

    orders = self._orders(orders)
    if 2 not in orders:
      raise ValueError('Cannot compute effective masses without second order term.')

    results = []
    breakpoints = self.breakpoints
    recipcell = inv(self.structure.cell).T * 2e0 * pi / self.structure.scale
    for start, end, direction in zip( breakpoints[:-1], 
                                      breakpoints[1:], 
                                      self.directions ):
      kpoints = self.kpoints[start:end]
      x = dot(direction, dot(recipcell, kpoints.T)) 
      measurements = self.eigenvalues[start:end].copy()
      parameters = concatenate([x[:, None]**i / factorial(i) for i in orders], axis=1)
      fit = lstsq(parameters, measurements)
      results.append(fit[0][orders.index(2)])

    result = (array(results) * self.eigenvalues.units * angstrom**2 / h_bar**2)
    return 1./result.rescale(1/emass)

  def fit_directions(self, orders=None):
    """ Returns fit for computed directions.
    
        When dealing with degenerate states, it is better to look at each
        computed direction separately, since the order of bands might depend on
        the direction (in which case it is difficult to construct a tensor).
    """
    from numpy import dot, concatenate, pi
    from numpy.linalg import inv, lstsq
    from math import factorial
    from ...error import ValueError

    orders = self._orders(orders)
    if 2 not in orders:
      raise ValueError('Cannot compute effective masses without second order term.')

    results = []
    breakpoints = self.breakpoints
    recipcell = inv(self.structure.cell).T * 2e0 * pi / self.structure.scale
    for start, end, direction in zip( breakpoints[:-1], 
                                      breakpoints[1:], 
                                      self.directions ):
      kpoints = self.kpoints[start:end]
      x = dot(direction, dot(recipcell, kpoints.T)) 
      measurements = self.eigenvalues[start:end].copy()
      parameters = concatenate([x[:, None]**i / factorial(i) for i in orders], axis=1)
      fit = lstsq(parameters, measurements)
      results.append(fit[0])

    return results

        

class _OnFinish(object):
  """ Called when effective mass calculation finishes. 
  
      Adds some data to the calculation so we can figure out what the arguments
      to the call.
  """
  def __init__(self, previous, outdir, center, nbpoints, range, directions):
    super(_OnFinish, self).__init__()
    self.previous = previous
    self.outdir = outdir
    self.center = [u for u in center]
    self.nbpoints = nbpoints
    self.directions = directions
    self.range = range
  def __call__(self, *args, **kwargs):
    from os.path import join
    # first calls previous onfinish.
    if self.previous is not None: self.previous(*args, **kwargs)
    # then adds data
    header = ''.join(['#']*20)
    with open(join(self.outdir, 'prop.out'), 'a') as file: 
      file.write('{0} {1} {0}\n'.format(header, 'EMASS DETAILS'))
      file.write('center={0!r}\n'.format(self.center))
      file.write('nbpoints={0!r}\n'.format(self.nbpoints))
      file.write('range={0!r}\n'.format(self.range))
      file.write('directions={0!r}\n'.format(self.directions))
      file.write('{0} END {1} {0}\n'.format(header, 'EMASS DETAILS'))

def iter_emass( functional, structure=None, outdir=None, workdir=None,
                center=None, nbpoints=5, range=0.1, directions=None, **kwargs ):
  """ Computes effective mass tensor for a given kpoint. 
  
      Computes Taylor expansion of a Hamiltonian around a given k-point.
      The expansion will be obtained by fitting an n-th order expansion to a
      set of eigenvalues in the neighberhood of the kpoint. In practice, a
      a fixed set of "band-structure" paths between
      diametrically opposed end-points are computed for a number of directions
      (001, 110, 111,...). The eigenvalues are obtained for each
      path using CRYSTAL_'s BAND keyword (in properties program). 

      :param functional:
          The effective masses for the given k-point are computed for a given
          fixed Hamiltonian. This fixed hamiltonian must first be generated
          from a self-consistent electronic minization calculation. For
          convenience, the fixed Hamiltonian can be either be given  on input
          or  computed within this function. As such, ``functional`` can be one
          of two types.

          - :py:class:`~pylada.dftcrystal.functional.Functional` object: a
             self-consistent calculation is first performed. The ``structure``
             argument cannot be None.
          - extraction object: e.g. as returned by a call to a
            :py:class:`~pylada.dftcrystal.functional.Functional` object. It
            should point to the self-consistent calculation on which the
            effective mass calculations will be based. The ``structure``
            argument is ignored.

      :param structure:
        If the argument ``functional`` is indeed a CRYSTAL functional, then
        this routine will first perform as self-consistent electronic
        minization (or whatever the functional does) in order to obtain a
        Hamiltonian. This argument is the structure for which to compute the
        Hamiltonian. If the argument ``functional`` is an extraction object,
        than ``structure`` is ignored.

      :param outdir: Output directory.
      :param workdir: Directory where calculations are performed.
      :param center:
          k-point for which to perform taylor expansion. Cartesian coordinates,
          units of 1/angstrom.
      :type  center: 3d-vector in fractional coordinates
      :param int nbpoints: 
          Number of points per "band-structure" lines.
      :param float range:
          Extent of the grid around the central k-point.
      :param directions:
          Array of directions (cartesian coordinates). If None, defaults to a
          reasonable set of directions: 001, 110, 111 and so forth. Note that
          if given on input, then the tensors should not be extracted. The
          directions are normalized. Eventually, the paths will extend from
          ``directions/norm(directions)*range`` to
          ``-directions/norm(directions)*range``.

      .. warning:: 
      
         In the case of degenerate states, it may be difficult to compute the
         effective mass tensors, especially in the case of warped eigenvalue
         surfaces, as in Si.
  """
  from numpy import array, zeros, sqrt, dot
  from numpy.linalg import norm
  from ...error import input as InputError
  from . import Properties

  # default kpoint is gamma.
  if center is None: center = zeros(3, dtype='float64')

  # If has an 'iter' function, then calls it. 
  if hasattr(functional, 'iter'): 
    if structure is None:
      raise InputError( 'If the first argument to iter_emass is a functional, '\
                        'then a structure must also be given on which to '     \
                        'apply the CRYSTAL functional.' )
    for input in functional.iter( structure, outdir=outdir, workdir=workdir,
                                  **kwargs ):
      if getattr(input, 'success', False): continue
      elif hasattr(input, 'success'):
        yield Extract(outdir)
        return
      yield input 
  # if is callable, then calls it.
  elif hasattr(functional, '__call__'):
    input = functional(structure, outdir=outdir, workdir=workdir, **kwargs)
  # otherwise, assume it is an extraction object.
  else: input = functional

  # check that self-consistent run was successful.
  if not input.success:
    yield Extract(outdir, input=input)
    return

  # creates Properties functional
  functional = Properties(input=input)
  # creates polar grid
  if directions is None:
    bands = array([ [1, 0, 0], [-1, 0, 0],
                    [0, 1, 0], [0, -1, 0],
                    [0, 0, 1], [0, 0, -1],
                    [1, 0, 1], [-1, 0, -1],
                    [0, 1, 1], [0, -1, -1],
                    [1, 1, 0], [-1, -1, 0],
                    [1, 0, -1], [-1, 0, 1],
                    [0, -1, 1], [0, 1, -1],
                    [-1, 1, 0], [1, -1, 0],
                    [1, 1, 1], [-1, -1, -1],
                    [1, 1, -1], [-1, -1, 1],
                    [1, -1, 1], [-1, 1, -1],
                    [-1, 1, 1], [1, -1, -1] ], dtype='float64')
    bands[6:18] *= 1e0/sqrt(2.)
    bands[18:] *= 1e0/sqrt(3.)
  else: 
    directions = array(directions).reshape(-1, 3)
    bands = array([d/norm(d) for d in directions])
    bands = array([(d, -d) for d in bands]).reshape(-1, 3)
  bands *= range
  functional.band = dot(input.structure.cell.T, (bands + center).T).T
  functional.band.nsub = len(functional.band) * nbpoints

  # now call functional.
  for u in functional.iter(outdir=outdir, workdir=workdir):
    if getattr(u, 'success', False): continue
    if hasattr(u, 'success'): yield u; return
    # modify onfinish so that call arguments are added to the output file.
    onfinish = _OnFinish( u.onfinish, outdir, center, nbpoints, range,
                          directions )
    u.onfinish = onfinish
    yield u
    
  yield Extract(outdir, input=input)

iter_emass.Extract = Extract
""" Extractor class for effective mass. """
EMass = makeclass( 'EMass', Properties, iter_emass, None,
                    module='pylada.dftcrystal.properties.emass',
                    doc='Functional form of the '                              \
                        ':py:class:`pylada.dftcrystal.properties.emass'          \
                        '.iter_emass` method.' )

# Function call to effective mass. No iterations. returns when calculations are
# done or fail.
effective_mass = makefunc( 'effective_mass', iter_emass,
                           'pylada.dftcrystal.properties.emass' )
