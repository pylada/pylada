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

""" Numerical energy derivatives. """
__docformat__ = "restructuredtext en"
from .kpoints import KPoints, _reduced_grids_factory

class DDPoints(KPoints):
  """ Points to compute Derivate of the Dispersion. """
  def __init__(self, direction, center = None, order = 2, nbpoints = 0, stepsize=1e-2, relax=True):
    """ Initializes the dispersion derivative k-point object. 


        :Parameters:
          direction : 3-tuple 
            Direction for which to compute derivatives.
          center : 3-tuple
            Point at which to take derivative.
          order : int
            Order of the derivative.
          nbpoints : int
            Number of points to use in computing derivative.
          stepsize : float
            Distance between interpolation points. Default = 1e-2.
            Units of ``2|pi|/a``, with ``a=structure.scale``.
          relax : bool
            Whether kpoint should be relaxed from input to output vff
            structure. True by default.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    super(DDPoints, self).__init__()
    self.direction = direction
    """ Direction for which to compute derivatives. """
    self.center = center
    """ Point at which to take derivative. """
    self.order = order
    """ Order of the derivative. """
    self.nbpoints = nbpoints
    """ Number of points to use in computing derivative. """
    self.stepsize = stepsize
    """ Distance between interpolation points. 

        Units of ``2|pi|/a``, with ``a=structure.scale``.
        
        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  @property
  def parameters(self): 
    """ Parameters for derivation. """
    assert hasattr(self, "_parameters"),\
           RuntimeError("parameters attribute cannot be accessed befor calculation.")
    return self._parameters

  def _mnk(self, input, output):
    """ Yields lines of k-points to perform numerical derivation. """
    from math import factorial
    from numpy import array, dot, zeros, pi
    from numpy.linalg import norm, inv
    from quantities import angstrom
    from ..physics import a0

    assert norm(self.direction) > 1e-12, ValueError("Direction cannot be null.")
    assert self.order > 0, ValueError("Derivation order cannot be zero.")
    assert self.stepsize > 0, ValueError("Stepsize must be positive.")
    assert abs(output.scale) > 1e-12, ValueError("scale is zero in structure.")
    assert self.order > 0, ValueError("Order of derivative should be positive.")

    nbpoints  = max(self.order+1, self.nbpoints)
    direction = array(self.direction, dtype="float64") 
    center    = array(self.center, dtype="float64") if self.center is not None\
                else zeros((3,), dtype="float64")
    if self.relax:
      deformation = dot(inv(output.cell.T), input.cell.T)
      direction = dot(deformation, direction)
      center = dot(deformation, center)
    direction /= norm(self.direction) # normalizes direction.

    # yields vector at point derivation for odd numbers of points.
    if nbpoints % 2 == 1: yield 1, center

    # yields all other points.
    start = 1 if nbpoints % 2 == 1 else 0.5
    parameters = zeros(shape=(nbpoints, self.order+1), dtype="float64") 
    parameters[:,0] = 1
    units = 2e0 * pi * a0.rescale(angstrom) / output.scale 
    for i in range(0, nbpoints - nbpoints%2): 
      # computes position on derivation line.
      s = (1 if i % 2 == 0 else -1) * self.stepsize * (i//2+start)
      # yields reciprocal space vector where to do calculation.
      yield 1, center + direction * s
      # sets up regression parameters.
      parameters[i + nbpoints%2, 1:] = [ pow(s*units, n)/float(factorial(n))\
                                         for n in range(1, self.order+1) ]

    # saves parameters for later use.
    self._parameters = parameters

  def __repr__(self):
    result = '{0.__class__.__name__}({1}'.format(self, repr(self.direction))
    do_key = self.center is None
    if self.center is not None: result += ', {0}'.format(repr(self.center))
    if self.order != 2:
      result += ', {1}{0}'.format(repr(self.relax), 'relax=' if do_key else '') 
    else: do_key = True
    if self.nbpoints != 0:
      result += ', {1}{0.nbpoints}'.format(self, 'nbpoints=' if do_key else '') 
    else: do_key = True
    if self.stepsize != 1e-2:
      result += ', {1}{0.stepsize}'.format(self, 'stepsize=' if do_key else '') 
    else: do_key = True
    return result + ')'

ReducedDDPoints    = _reduced_grids_factory('ReducedDDPoints', DDPoints)


class ChainedDDPoints(KPoints):
  """ Chains together different directions.
  
      The points is to only use those calculations which are necessary.
  """
  def __init__(self, direction = None, *args, **kwargs):
    """ Initializes a set of chained directions. 
    
        The first argument, ``direction``, can now be a list of directions. 
    """
    from numpy import array
    super(ChainedDDPoints, self).__init__()
    direction = array(direction, dtype="float64")
    if direction.ndim == 1: direction = array([direction])
    self.ddpoints = [DDPoints(dir, *args, **kwargs) for dir in direction]
    """ List of DDPoints instances to chain. """

  @property
  def center(self):
    """ Point at which to take derivative. """
    return self.ddpoints[0].center
  @center.setter
  def center(self, value):
    from numpy import array
    value = array(value, dtype='float64')
    for u in self.ddpoints: u.center = value

  @property
  def order(self):
    """ Order of the derivative. """
    return self.ddpoints[0].order
  @order.setter
  def order(self, value):
    from numpy import array
    value = array(value, dtype='float64')
    for u in self.ddpoints: u.order = value

  @property
  def nbpoints(self):
    """ Number of points to use in computing derivative. """
    return self.ddpoints[0].nbpoints
  @nbpoints.setter
  def nbpoints(self, value):
    for u in self.ddpoints: u.nbpoints = value

  @property
  def stepsize(self):
    """ Distance between interpolation points. 

        Units of ``2|pi|/a``, with ``a=structure.scale``.
        
        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    return self.ddpoints[0].stepsize
  @stepsize.setter
  def stepsize(self, value):
    for u in self.ddpoints: u.stepsize = value

  @property
  def relax(self):
    """ Whether to deform kpoints to the relaxed structure. """
    return self.ddpoints[0].relax
  @relax.setter
  def relax(self, value):
    for u in self.ddpoints: u.stepsize = value

  @property
  def directions(self):
    """ List of all directions. """
    from numpy import array
    return array([u.direction for u in self.ddpoints])

  def __repr__(self): 
    """ Returns string representing this object. """
    from copy import copy
    a = copy(self.ddpoints[0])
    a.direction = self.directions 
    return repr(self.ddpoints[0]).replace("DDPoints", self.__class__.__name__)

  @property
  def parameters(self):
    """ List of parameters for each direction. """
    from numpy import array
    return array([u.parameters for u in self.ddpoints])

  def _mnk(self, input, output):
    """ Yields lines of k-points to perform numerical derivation. """
    for point in self.ddpoints:
      for result in point._mnk(input, output): 
        yield result


ReducedChainedDDPoints  = _reduced_grids_factory('ReducedChainedDDPoints', ChainedDDPoints)



def reciprocal( escan, structure, outdir = None, comm = None, direction=(0,0,1), order = 1, \
                nbpoints = None, stepsize = 1e-2, center = None, lstsq = None, **kwargs ):
  """ Computes effective mass for a given direction.

      :Parameters:
        escan : `Escan` or `KEscan`
          Emiprical pseudo-potential functional.
        structure : `crystal.Structure`
          The structure for wich to compute effective masses.
        outdir : str
          Directory where to save results of calculation.
        comm : `pylada.mpi.Communicator` or None
          MPI communicator containing processes with which to perform
          calculation.
        direction : 3-tuple or list of 3-tuple
          direction for which to compute derivatives.
        order : int
          Highest order derivative to perform. Defaults to 1.
        nbpoints : int
          Number of points with wich to compute derivatives.
          Should be at least order + 1. Default = order + 1. 
        stepsize : float
          Distance between interpolation points. Default = 1e-2.
          Units of ``2|pi|/a``, with ``a=structure.scale``.
        center : 3-tuple
          k-point where to take derivative. Units of ``2|pi|/a``, with
          ``a=structure.scale``.
        lstsq 
          Linear least square method. The first two parameters should
          be same as numpy.linalg.lstsq. Other parameters can be passed as extra
          parameters. Defaults to numpy.linalg.lstsq.
        kwargs 
          Extra parameters which are passed on first to escan (if escan
          object has an attribute with the same name), then to the linear least
          square fit method. Note that this will *not* change the external escan
          object.  This function is stateless. 

      :return: Same as return from lstsq.

      :warning: escan.nbstates (or nbstates if passed as arg) must be the
                *spin-polarized* number of electrons, whatever escan thinks.

      .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI

  """
  from numpy import array, sort
  from numpy.linalg import lstsq as np_lstsq
  from quantities import hartree
# from .kescan import KEscan

  # takes care of default parameters.
# if not isinstance(escan, KEscan): escan = KEscan(escan=escan)
  if center is None: center = kwargs.pop("kpoint", escan.kpoint)
  center = array(center, dtype="float64")
  relax = kwargs.pop("do_relax_kpoint", escan.do_relax_kpoint)
  if outdir is None: outdir = "reciprocal"
  if lstsq is None: lstsq = np_lstsq
  direction = array(direction, dtype="float64")

  # creates kpoints object.
  if direction.ndim == 2:
    kpoints = ReducedChainedDDPoints(direction, center, order, nbpoints, stepsize, relax)
  else: 
    kpoints = ReducedDDPoints(direction, center, order, nbpoints, stepsize, relax)

  # performs calculations.
  out = escan(structure, outdir=outdir, comm=comm, kpoints=kpoints, **kwargs)
  # makes sure we have parameters. 
  for k in kpoints(out.input_structure, out.structure): continue

  # sorts eigenvalues at each kpoint and rescales to hartree.
  measurements = sort(out.eigenvalues.rescale(hartree), axis=1) 
  
  # finally, performs least-square fit and returns everything.
  if direction.ndim == 2: # for each direction.
    result = []
    nbpoints = max(kpoints.order+1, kpoints.nbpoints)
    for i in xrange(direction.shape[0]):
      j, k = i * nbpoints, (i+1) * nbpoints
      result.append(lstsq(kpoints.parameters[i,:,:], measurements[j:k, :]))
    return result
  # or for the lone direction.
  return lstsq( kpoints.parameters, measurements )

