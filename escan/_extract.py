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

""" Module to extract esca and vff ouput. """
__docformat__ = "restructuredtext en"
__all__ = ['Extract']

from ..opt.decorators import broadcast_result, make_cached, FileCache
from ..opt import AbstractExtractBase, OutcarSearchMixin

class Extract(AbstractExtractBase, OutcarSearchMixin):
  """ A class to extract data from ESCAN output files. 
  
      This class helps to extract information from the escan output, including
      escan and vff parameters, relaxed crystal structure, eigenvalues, and
      optionally, wavefunctions in real and reciprocal space. Where possible
      quantities have attached units, using the package ``quantities``.
  """
  def __init__(self, directory = None, comm = None, escan = None):
    """ Initializes ESCAN extraction class. 
    
        :Parameters: 
          directory : str or None
            Directory where escan output files are located. Defaults to current
            working directory if None.
          comm : `mpi.Communicator`
            Communicator containing as many processes as were used to perform
            calculations. This is only mandatory when using wavefunctions in
            some way.
          escan : pylada.escan.Escan
            Wrapper around the escan functional.
    """
    from ..vff import Extract as VffExtract
    from . import Escan

    super(Extract, self).__init__(directory=directory, comm=None)

    if escan is None: escan = Escan()
    
    self.OUTCAR = escan.OUTCAR
    """ OUTCAR file to extract stuff from. """
    self.FUNCCAR = escan._FUNCCAR
    """ Pickle to FUNCCAR. """
    
    self._vffout = VffExtract(directory, comm = comm, vff = escan.vff)
    """ Private reference to vff extraction object. """

  def __funccar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.FUNCCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  @property 
  def comm(self):
    """ Communicator over which to sync output. """
    return self._vffout.comm 
  @comm.setter
  def comm(self, value):
    from ..mpi import Communicator
    if hasattr(self, '_vffout'): self._vffout.comm = Communicator(value)

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def kpoint(self):
    """ K-point in this calculation. """
    from numpy import array, zeros, pi
    from quantities import angstrom
    from ..physics import a0
    from re import M
    result = self._find_first_OUTCAR(r'\s*ikpt,akx,aky,akz\s+(\d+)' + 3*r'\s+(\S+)', M)
    assert result is not None,\
           RuntimeError('Could not find kpoint in file {0};'.format(self.__outcar__().name))
    if result.group(1) == '0': return zeros((3,), dtype='float64')
    return array([result.group(2), result.group(3), result.group(4)], dtype='float64') / 2e0 / pi\
           * (self.solo().structure.scale * angstrom).rescale(a0).magnitude

  
  def __directory_hook__(self):
    """ Called whenever the directory changes. """
    super(Extract, self).__directory_hook__()
    self._vffout.directory = self.directory

  def uncache(self): 
    """ Uncache values. """
    super(Extract, self).uncache(self)
    self._vffout.uncache()

  def __copy__(self):
    """ Returns a shallow copy of this object. """
    result = super(Extract, self).__copy__()
    result._vffout = self._vffout.__copy__()
    return result

  def copy(self, **kwargs):
    """ Returns a shallow copy of this object.

        :param kwargs:
          Any keyword argument is set as an attribute of this object.
    """
    result = self.__copy__()
    if 'comm' in kwargs: result.comm = kwargs.pop('comm')
    for k, v in kwargs.items():
      if hasattr(result._vffout, k):
        setattr(result._vffout, k, v)
        if hasattr(result, k): setattr(result, k, v)
      else: setattr(result, k, v)
    return result

  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for Escan success.
        
        At this point, checks for files and 
    """
    from numpy import any, isnan
    from re import compile
    do_escan_re = compile(r'functional\.do_escan\s*=')

    try:
      good = 0
      is_do_escan = True
      with self.__outcar__() as file:
        for line in file:
          if line.find("FINAL eigen energies, in eV") != -1: good += 1
          elif do_escan_re.search(line) is not None: is_do_escan = eval(line.split()[-1])
          elif line.find("# Computed ESCAN in:") != -1: good += 1; break
      if good == 1 and not is_do_escan: return True
      if good == 2 and is_do_escan: return not any(isnan(self.solo().eigenvalues))
      return False
    except: return False

  @property
  @make_cached
  def functional(self):
    """ Greps escan functional from OUTCAR. """
    from cPickle import load
    from ..vff import _get_script_text
    from . import exec_input
    
    # tries to read from pickle.
    try:
      with self.__funccar__() as file: result = load(file)
    except: pass 
    else: return result


    # tries to read from outcar.
    @broadcast_result(attr=True, which=0)
    def get_functional(this):
      with self.__outcar__() as file: return _get_script_text(file, "Escan")
    # moves to output directory to get relative paths right.
    input = exec_input(get_functional(self), {'vff_functional': self.vff})
    return input.escan_functional if "escan_functional" in input.__dict__\
           else input.functional

  @property 
  def escan(self): 
    """ Alias for functional. """
    from warnings import warn
    warn(DeprecationWarning('escan attribute is deprecated in favor of functional.'), stacklevel=2)
    return self.functional

  @property
  def _double_trouble(self):
    """ Returns true, if non-spin polarized or Kammer calculations. """
    from numpy.linalg import norm
    from . import soH
    seul = self.solo()
    if seul.functional.nbstates  ==   1: return False
    if seul.functional.potential != soH: return True
    return norm(seul.functional.kpoint) < 1e-12


  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def eigenvalues(self):
    """ Greps eigenvalues from OUTCAR. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    from quantities import eV
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file {0}.".format(path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL eigen energies, in eV") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find("*********************************") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    if self._double_trouble: result = [result[i/2] for i in range(2*len(result))]
    return array(result, dtype="float64") * eV

  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def convergence(self):
    """ Greps eigenvalue convergence errors from OUTCAR. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL err of each states, A.U") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find(" FINAL eigen energies, in eV") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    if self._double_trouble: result = [result[i/2] for i in range(2*len(result))]
    return array(result, dtype="float64") 

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nnodes(self):
    """ Greps eigenvalue convergence errors from OUTCAR. """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" nnodes =") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
    return int(line.split()[2])


  @property
  @make_cached
  def gwfns(self):
    """ Creates list of Wavefuntion objects. """
    from _wfns import Wavefunction
    result = []
    if self.is_spinor:
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i/2,0],\
                                        self.raw_gwfns[:,i/2,1], attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(self.comm, i, eig,\
                                        -self.raw_gwfns[self.inverse_indices,i/2,1].conjugate(),\
                                         self.raw_gwfns[self.inverse_indices,i/2,0].conjugate(), \
                                        attenuation = self.attenuation) )
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i,0],\
                                      self.raw_gwfns[:,i,1], attenuation = self.attenuation) )
    else: # no spin polarization.
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i/2,0],\
                                        attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(self.comm, i, eig, \
                                        self.raw_gwfns[self.inverse_indices,i/2,0], \
                                        attenuation = self.attenuation) )
          result.append(result[-1])
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i,0],\
                                      None, self.attenuation) )
          result.append(result[-1])
    return result

  @property
  @make_cached
  def rwfns(self):
    """ Creates list of rWavefuntion objects. """
    from ._wfns import rWavefunction, gtor_fourrier
    result = []
    if self.is_spinor:
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          if i%2 == 0:
            rwfn = rWavefunction( self.comm, i, eig, self._raw_rwfns[:,i/2,0],\
                                  self._raw_rwfns[:,i/2,1])
          else: 
            rwfn = rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i/2,1].conjugate(),\
                                 -self._raw_rwfns[:,i/2,0].conjugate())
          result.append(rwfn)
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          rwfn = rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i,0], self._raw_rwfns[:,i,1])
          result.append(rwfn)
    else: # no spin polarization.
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          if i%2 == 0: 
            result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i/2,0]) )
          else:
            result.append( rWavefunction( self.comm, i, eig,
                                          -self._raw_rwfns[:,i/2,0].conjugate()) )
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i,0]) )
    return result

  @property
  def fft_mesh(self):
    """ Returns gvector mesh. """
    regex = self._find_first_OUTCAR(r"n1,n2,n3\s*=\s*(\d+)\s*(\d+)\s*(\d+)")
    return int(regex.group(1)), int(regex.group(2)), int(regex.group(3))

  @property
  def raw_rwfns(self):
    """ Raw real-space wavefunction data. """
    if not hasattr(self, "_raw_rwfns"): self.rwfns # creates data
    return self._raw_rwfns

  @property
  @make_cached
  def _raw_gwfns_data(self):
    """ Reads and caches g-space wavefunction data. 
    
        This property is a tuple holding information about the wavefunctions.
        
        - a spin by N by x matrix holding the N wavefuntions/spinor.
        - a 3 by x matrix with each row a G-vector in units of
          `pylada.physics.reduced_reciprocal_au`.
        - a 3 by x matrix with each row a R-vector in atomic units.
        - one-dimensional array of real coefficients to smooth higher energy G-vectors.
        - one-dimensional array of integer indices to map G-vectors to -G.
    """
    # check for mpi first
    from .. import pylada_with_mpi
    assert pylada_with_mpi, RuntimeError("Pylada loaded without mpi. Cannot read wavefunctions.")
    # then check for function.
    from os.path import exists
    from numpy import sqrt, abs, array
    from numpy.linalg import norm, det
    from quantities import angstrom, pi
    from ..opt import redirect
    from ..opt.changedir import Changedir
    from ..physics import a0, reduced_reciprocal_au
    from ._escan import read_wavefunctions
    from . import soH

    assert self.comm.real, ValueError("MPI needed to play with wavefunctions.")
    assert self.success, RuntimeError("Run was unsuccessful. Cannot read wavefunctions.")
    assert self.comm.size >= self.nnodes or norm(self.kpoint) < 1e-12,\
           RuntimeError( "Cannot read wavefunctions with fewer procs "\
                         "({0} < {1}) than written when not at Gamma."\
                         .format(self.comm.size, self.nnodes) )

    # case where we need same number  of procs for reading as writing.
    if norm(self.kpoint) > 1e-12 and self.comm.size != self.nnodes:
      local_comm = self.comm.split(self.comm.rank < self.nnodes)
      return self.copy(comm=local_comm)._raw_gwfns_data if self.comm.rank < self.nnodes \
             else (None, None, None, None, None)
                          
    kpoint = array(self.functional.kpoint, dtype="float64")
    scale = self.structure.scale / float(a0.rescale(angstrom))
    with Changedir(self.directory, comm=self.comm) as directory:
      assert exists(self.functional.WAVECAR),\
             IOError("{0} does not exist.".format(self.functional.WAVECAR))
      nbstates = self.functional.nbstates
      if self.functional.potential != soH or norm(self.functional.kpoint) < 1e-6:
        nbstates = max(1, self.nbstates/2)
      with redirect(fout="") as streams:
        result = read_wavefunctions(self.functional, range(nbstates), kpoint, scale, self.comm)
    self.comm.barrier()

    cell = self.structure.cell * self.structure.scale * angstrom
    normalization = abs(det(cell.rescale(a0)))
    return result[0] * sqrt(normalization), result[1] * 0.5 / pi * reduced_reciprocal_au,\
           result[2] * a0, result[3], result[4]

  @property
  def raw_gwfns(self):
    """ Raw wavefunction data in g-space. 
    
        Numpy array with three axis: (i) g-vectors, (ii) bands, (iii) spins:

        >>> self.raw_gwfns[:,0, 0] # spin up components of band index 0.
        >>> self.raw_gwfns[:,0, 1] # spin down components of band index 0.
        >>> for i, g in enumerate(self.gvectors): # looks for G=0 component
        >>>   if np.linalg.norm(g) < 1e-8: break
        >>> self.raw_gwfns[i, 0, 0] # G=0 component of spin-up wavefunction with band-index 0.

        The band index is the one from ESCAN, eg. it is different if user
        Krammer doubling or not, etc. This data is exactly as read from disk.
    """
    return self._raw_gwfns_data[0]

  @property
  def gvectors(self):
    """ G-vector values of wavefuntions. """
    return self._raw_gwfns_data[1]

  @property
  def rvectors(self):
    """ R-vector values of wavefuntions. """
    return self._raw_gwfns_data[2]

  @property
  def attenuation(self):
    """ G-vector attenuation values of wavefuntions. """
    return self._raw_gwfns_data[3]

  @property
  def inverse_indices(self):
    """ Indices to -G vectors of wavefuntions. """
    return self._raw_gwfns_data[4]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _wavefunction_path(self): return self.solo().functional.WAVECAR

  @property
  def is_spinor(self):
    """ True if wavefunction is a spinor. """
    from . import soH
    return self.functional.potential == soH

  @property
  def is_krammer(self):
    """ True if wavefunction has krammer degenerate equivalent. """
    return self.functional.is_krammer

  @property
  def vff(self):
    """ Vff functional. """
    return self._vffout.functional

  @FileCache('DIPOLESCAR')
  def _dipoles(self, attenuate=False):
    """ Computes dipole matrix element between all states.
    
        This routine caches results in a file. The routine above should check
        that the arguments are the same.
    """
    from numpy import zeros
    result = zeros(shape=(len(self.eigenvalues), len(self.eigenvalues), 3), dtype="complex64")
    gvectors = self.gvectors
    for i, wfnA in enumerate(self.gwfns):
      for j, wfnB in enumerate(self.gwfns):
        if j <= i: continue
        result[i, j, :] = wfnA.braket(gvectors, wfnB, attenuate=attenuate)
    for i in range(len(self.gwfns)):
      for j in range(len(self.gwfns)):
        if j < i: result[i, j, :] = result[j, i, :].conjugate()
    result = result * self.gvectors.units
    return attenuate, result


  @FileCache('EMASSCAR')
  def effective_mass_tensor(self, attenuate=False, degeneracy=1e-8):
    """ Returns list of effective mass tensor (1/m_xy) for each band. """
    from numpy import zeros, outer, identity
    from pylada.physics import electronic_mass, h_bar

    # gets dipoles.
    dipoles = zeros((len(self.eigenvalues), len(self.eigenvalues), 3), dtype="complex64") #self.dipoles(attenuate)
    gvectors = self.gvectors
    for i, wfnA in enumerate(self.gwfns):
      for j, wfnB in enumerate(self.gwfns):
        dipoles[i, j, :] = wfnA.braket(gvectors, wfnB, attenuate=attenuate)
    dipoles = dipoles * gvectors.units

    # then computes second order terms.
    result = zeros(shape=(len(self.eigenvalues), 3,3), dtype="float64")
    for i, eigA in enumerate(self.eigenvalues):
      for j, eigB in enumerate(self.eigenvalues):
        if abs(eigA - eigB) < degeneracy: continue
        vector = (dipoles[j, i, :] - dipoles[i, i, :].real) * h_bar**2/electronic_mass
        result[i, :, :] += 2 * outer(vector.conjugate(), vector).real / (eigA - eigB)
    
    units = (dipoles.units * h_bar**2 / electronic_mass)**2 / self.eigenvalues.units / h_bar**2
    return identity(3)/electronic_mass + result * units

  def dipoles(self, attenuate=False):
    """ Computes dipole matrix element between vbm and cbm. """
    # gets result, possibly from cache file.
    a2, result = self._dipoles(attenuate)

    uncache = a2 != attenuate
    if uncache: 
      from os.path import join
      from os import remove
      try: remove(join(self.directory, "DIPOLESCAR"))
      except: pass
      return self._dipoles(attenuate)[-1]
    return result

# def effective_mass_tensor(self, attenuate=False, degeneracy=1e-3):
#   """ Computes effective masses tensor using to dipole matrix elements. """
#   from numpy import array, dot, multiply, sum, mean, zeros, identity
#   from numpy.linalg import inv
#   from quantities import dimensionless
#   from ..physics import electronic_mass, h_bar

#   # compute dipoles.
#   dipoles = self.dipoles(attenuate)
#   # create equivalent array with units/(e0 - e1)
#   units = 2e0 * h_bar**2 / electronic_mass 
#   factor = array([ [ (units/(e0 -e1) if abs(e0 - e1) > 1e-8 else 0) for e1 in self.eigenvalues ]\
#                    for e0 in self.eigenvalues ]) * units.units / self.eigenvalues.units

#   # In case of degenerate subspaces, we should average of different values. 
#   # First figure out what the degeneracy classes are.
#   classes = [[0]]
#   for i, (first, second) in enumerate(zip(self.eigenvalues[:-1], self.eigenvalues[1:])): 
#     if abs(first - second) < degeneracy and degeneracy > 0e0: classes[-1].append(i+1)
#     else: classes.append([i+1])
#   classes = [r for r in classes if len(r) > 1]

#   # now computes all oscillator strengths.
#   result = zeros((3,3,len(self.eigenvalues)), dtype="float64")
#   for i in xrange(3):
#     x = zeros((3,),dtype="float64"); x[i] = 1e0
#     px = dot(dipoles, x)
#     for j in xrange(3):
#       y = zeros((3,),dtype="float64"); y[j] = 1e0
#       py = dot(dipoles, y)
#       dummy = multiply(multiply(px.conjugate(), py).real, factor).simplified 
#       assert dummy.units == dimensionless
#       dummy = sum(dummy, axis=1)
#       for class_ in classes:
#         dummy[class_] = mean(dummy[class_])
#         dummy[class_] = mean(dummy[class_])
#       result[i,j,:] = dummy

#   # and return result.
#   for i in range(len(self.eigenvalues)):
#     result[:,:, i] = 1/(result[:, :, i] + identity(3))
#   return electronic_mass * result

  def __getattr__(self, name):
    """ Passes on public attributes to vff extractor, then to escan functional. """
    if name[0] != '_':
      if name in dir(self._vffout): return getattr(self._vffout, name)
      elif self.success and hasattr(self.functional, name):
        return getattr(self.functional, name)
    raise AttributeError("Unknown attribute {0}. It could be the "\
                         "run is unsuccessfull.".format(name))

  def __dir__(self):
    """ Returns list attributes.
    
        Since __getattr__ is modified, we need to make sure __dir__ returns a
        complete list of attributes. This is usefull  for command-line
        completion in ipython.
    """
    exclude = set(["add_potential", "write_escan_input"])
    result = [u for u in self.__dict__.keys() if u[0] != "_"]
    result.extend( [u for u in dir(self.__class__) if u[0] != "_"] )
    result.extend( [u for u in dir(self._vffout) if u[0] != "_"] )
    if self.success: result.extend( [u for u in dir(self.functional) if u[0] != "_"] )
    return list( set(result) - exclude )

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        :kwarg errors: Include stderr files.
        :type errors: bool
        :kwarg wavecar: Include wavefunctions.
        :type wavecar: bool
        :kwarg potcar: Include potential.
        :type potcar: bool
        :kwarg poscar: Include atom.config file.
        :type poscar: bool
    """
    from os.path import join, exists
    files = [self.OUTCAR, self.FUNCCAR, "DIPOLESCAR", "EMASSCAR"]
    if kwargs.get("poscar", False):
      try: files.append(self.functional._POSCAR)
      except: pass
    if kwargs.get("wavecar", False):
      try: files.append(self.functional.WAVECAR)
      except: pass
    if kwargs.get("errors", False):
      try: files.append(self.functional.ERRCAR)
      except: pass
    if kwargs.get("potcar", False):
      try: files.append(self.functional._POTCAR)
      except: pass
    for file in files:
      file = join(self.directory, file)
      if exists(file): yield file
    for file in self._vffout.iterfiles(**kwargs): yield file
    
