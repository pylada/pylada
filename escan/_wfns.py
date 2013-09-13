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

""" Classes to handle ESCAN wavefunctions. """
__docformat__  = 'restructuredtext en'

class rWavefunction(object):
  """ Defines the wavefunction of a single band in realspace. """
  is_gspace = False
  """ False since this is a r-space wavefunction. """

  def __init__(self, comm, index, eigenvalue, up, down = None):
    """ Initializes a spinor. """
    self.comm = comm
    """ Communicator over which wavefunctions are spread. """
    self.index = index
    """ Escan index of the wavefunction. """
    self.eigenvalue = eigenvalue
    """ Eigenvalue corresponding to this eigenvector. """
    self.up = up
    """ Up component (or unpolarized component) of the spinor. """
    self.down = down
    """ Down component of the spinor. None if not a spinor. """

  def expectation_value(self, operator):
    """ Returns expectation value of this operator. 
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    return self.braket(operator, self)

  def braket(self, operator, ket):
    """ Returns <this wavefunction|operator|ket wavefuntion>.
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply, transpose
    a = conjugate(self.up)
    if operator is None: b = ket.up
    elif operator.shape[0] == ket.up.shape[0] and operator.shape[-1] != ket.up.shape[0]:
      b = multiply(transpose(operator), ket.up) 
    else: b = multiply(operator, ket.up) 
    result = dot(b, a)
    if self.down is not None: 
      a = conjugate(self.down)
      if operator is None: b = ket.down
      elif operator.shape[0] == ket.down.shape[0] and operator.shape[-1] != ket.down.shape[0]:
        b = multiply(transpose(operator), ket.down) 
      else: b = multiply(operator, ket.down) 
      result += dot(b, a)
    if self.comm.is_mpi: result = self.comm.all_reduce(result, lambda x,y: x+y)
    return result * operator.units if hasattr(operator, "units") else result

  @property
  def density(self):
    """ Returns partial density. """
    from numpy import multiply, conjugate
    result = multiply(self.up, conjugate(self.up))
    if self.down is not None: result += multiply(self.down, conjugate(self.down))
    return result.real


class Wavefunction(rWavefunction):
  is_gspace = True
  """ True since this is a g-space wavefunction. """

  def __init__(self, *args, **kwargs):
    """ Initializes a spinor. """
    self.attenuation = kwargs["attenuation"] if "attenuation" in kwargs else None
    """ Attenuation coefficients of high-energy G-vectors. """
    if "attenuation" in kwargs: del kwargs["attenuation"]
    super(Wavefunction, self).__init__(*args, **kwargs)

  def expectation_value(self, operator, attenuate = True):
    """ Returns expectation value of this operator. 
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    return self.braket(operator, self, attenuate=attenuate)

  def braket(self, operator, ket, attenuate = True):
    """ Returns <this wavefunction|operator|ket wavefuntion>.
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import multiply, transpose
    if attenuate is None:  a = operator
    elif operator is None: a = multiply(self.attenuation, self.attenuation)
    elif     operator.shape[0] == self.attenuation.shape[0] \
         and operator.shape[-1] != self.attenuation.shape[0]:
      a = multiply(transpose(operator), multiply(self.attenuation, self.attenuation))
    else: a = multiply(operator, multiply(self.attenuation, self.attenuation))
    return super(Wavefunction, self).braket(a, ket)


def gtor_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  """ g-space to r-space fourrier transform of wavefunctions.
  
      :Parameters:
        wavefunctions : numpy array
          Arry of wavefunctions.
        rvectors : numpy array
          Two-dimensional array of r-space vectors, with each row a position.
          The return r-space wavefunctions will be given with respect to these
          points. Each process should have different r-space vectors. Otherwise
          use another implementation.
        gvectors : numpy array
          Two-dimensional array of g-space vectors, with each row a
          (g-)position. The input wavefunctions should be given with respect to
          these points, in the same order, etc.
        comm : `mpi.Communicator`
          communicator over which the wavefunctions are distributed.  The
          return wavefunctions will also be dirstributed over these processes.
        axis : integer
          axis over which the wavefunctions are deployed, eg axis of
          `wavefunctions` which corresponds to `gvectors`. Independent
          (though simultaneous, implementation wise) fourrier transform will be
          performed over this axis for all other axis. 0 by default.
  """
  import numpy as np

  # serial version
  if not comm.is_mpi:
    v = np.tensordot(rvectors, gvectors, ((1),(1))) 
    v = np.exp(-1j * v * (rvectors.units * gvectors.units).simplified)
    return np.tensordot(v, wavefunctions, ((1),(axis))) / np.sqrt(float(len(rvectors)))

  # mpi version
  else: 
    result = None
    for node in range(comm.size):
      # sends rvectors from node to all
      r = comm.broadcast(rvectors, node)
      # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
      v = np.tensordot(r, gvectors, ((1),(1)))
      v = np.exp(-1j * v * (rvectors.units * gvectors.units).simplified)
      # computes fourrier transform for all wavefunctions simultaneously.
      dummy = np.tensordot(v, wavefunctions, ((1),(axis)))
      # reduce across processes
      if node == comm.rank: result = comm.reduce(dummy, lambda x,y: x+y, node)
      else: comm.reduce(dummy, lambda x,y: x+y, node)

    assert not np.any(np.isnan(result))
    norm = comm.all_reduce(len(rvectors), lambda x, y: x+y)
    return result / np.sqrt(float(norm))


def rtog_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  """ r-space to g-space fourrier transform of wavefunctions.
  
      :Parameters:
        wavefunctions : numpy array
          Arry of wavefunctions.
        rvectors : numpy array
          Two-dimensional array of r-space vectors, with each row a position.
          The return r-space wavefunctions will be given with respect to these
          points. Each process should have different r-space vectors. Otherwise
          use another implementation.
        gvectors : numpy array
          Two-dimensional array of g-space vectors, with each row a
          (g-)position. The input wavefunctions should be given with respect to
          these points, in the same order, etc.
        comm : `mpi.Communicator`
          communicator over which the wavefunctions are distributed.  The
          return wavefunctions will also be dirstributed over these processes.
        axis : integer
          axis over which the wavefunctions are deployed, eg axis of
          `wavefunctions` which corresponds to `gvectors`. Independent
          (though simultaneous, implementation wise) fourrier transform will be
          performed over this axis for all other axis. 0 by default.
  """
  import numpy as np

  assert not np.any(np.isnan(wavefunctions))
  assert not np.any(np.isnan(gvectors))
  assert not np.any(np.isnan(rvectors))
  result = None
  
  if not comm.is_mpi:
    v = np.tensordot(rvectors, gvectors, ((1),(1)))
    v = np.exp(1j * v * (rvectors.units * gvectors.units).simplified)
    last = wavefunctions.ndim-1
    result = np.dot(wavefunctions.swapaxes(axis, last), v).swapaxes(last, axis)
    return result/ float(result.shape[0])

  # mpi version
  else: 
    for node in range(comm.size):
      # sends rvectors from node to all
      g = comm.broadcast(gvectors, node)
      # computes all exponentials exp(-i r.g), with g in first dim, and r in second.
      v = np.tensordot(rvectors, g, ((1),(1)))
      v = np.exp(1j * v * (rvectors.units * gvectors.units).simplified)
      # computes fourrier transform for all wavefunctions simultaneously.
      # somehow, there is a problem with tensordot leading to nan numbers...
      assert not np.any(np.isnan(v))
      last = wavefunctions.ndim-1
      dummy = np.dot(wavefunctions.swapaxes(axis, last), v).swapaxes(last, axis)
      # reduce across processes
      assert not np.any(np.isnan(dummy))
      if node == comm.rank: result = comm.reduce(dummy, lambda x,y: x+y, node)
      else: comm.reduce(dummy, lambda x,y: x+y, node)


    assert not np.any(np.isnan(result))
    # gets normalization factor.
    norm = comm.all_reduce(rvectors.shape[0], lambda x,y:x+y)
    assert norm != 0
    return result/ float(norm)
