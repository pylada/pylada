""" Module to extract esca and vff ouput. """
__docformat__ = "restructuredtext en"

from ..vff import Extract as VffExtract, _get_script_text
from ..opt.decorators import broadcast_result, make_cached

class Extract(object):
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
          comm : boost.mpi.communicator
            Communicator containing as many processes as were used to perform
            calculations. This is only mandatory when using wavefunctions in
            some way.
          escan : lada.escan.Escan
            Wrapper around the escan functional.
    """
    from os import getcwd
    from . import Escan

    if escan == None: escan = Escan()
    super(Extract, self).__init__()
    
    self._vffout = VffExtract(directory, comm = comm, vff = escan.vff)
    """ Private reference to vff extraction object. """
    self.OUTCAR = escan.OUTCAR
    """ OUTCAR file to extract stuff from. """
    self.FUNCCAR = escan._FUNCCAR
    """ Pickle to FUNCCAR. """
    self.comm = comm
    """ Communicator for extracting stuff. 

        All procs will get same results at end of extraction. 
        Program will hang if not all procs are called when extracting some
        value. Instead, use `solo`.

        >>> extract.success # Ok
        >>> if comm.rank == 0: extract.success # will hang if comm.size != 1
        >>> if comm.rank == 0: extract.solo().success # Ok
    """
    self.directory = directory if directory != None else getcwd()
  
  def _get_directory(self):
    """ Directory with VASP output files """
    return self._directory
  def _set_directory(self, dir):
    from os.path import abspath, expanduser
    dir = abspath(expanduser(dir))
    if hasattr(self, "_directory"): 
      if dir != self._directory: self.uncache()
    self._directory = dir
  directory = property(_get_directory, _set_directory)

  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for Escan success.
        
        At this point, checks for files and 
    """
    from os.path import exists, join
    path = join(self.directory, self.OUTCAR) if len(self.directory) else self.OUTCAR
    if not exists(path): return False

    good = 0
    is_do_escan = True
    with open(path, "r") as file:
      for line in file:
        if line.find("FINAL eigen energies, in eV") != -1: good += 1
        if line.find("functional.do_escan              =") != -1:
          is_do_escan = eval(line.split()[-1])
        if line.find("# Computed ESCAN in:") != -1: good += 1; break
    return (good == 2 and is_do_escan) or (good == 1 and not is_do_escan)

  @property
  @make_cached
  def escan(self):
    """ Greps escan functional from OUTCAR. """
    from os.path import exists, join
    from numpy import array
    from cPickle import load
    from ..opt.changedir import Changedir
    from . import Escan, localH, nonlocalH, soH, AtomicPotential
    
    # tries to read from pickle.
    path = self.FUNCCAR
    if len(self.directory): path = join(self.directory, self.FUNCCAR)
    if exists(path):
      try:
        with open(path, "r") as file: result = load(file)
      except: pass 
      else: return result

    # tries to read from outcar.
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))

    @broadcast_result(attr=True, which=0)
    def get_functional(this):
      with open(path, "r") as file: return _get_script_text(file, "Escan")
    local_dict = { "lattice": self.lattice, "minimizer": self.minimizer,\
                   "vff_functional": self.vff, "Escan": Escan, "localH": localH,\
                   "nonlocalH": nonlocalH, "soH": soH, "AtomicPotential":AtomicPotential,\
                   "array": array }
    # moves to output directory to get relative paths right.
    with Changedir(self.directory, comm=self.comm) as cwd:
      exec get_functional(self) in globals(), local_dict
    return local_dict["escan_functional"] if "escan_functional" in local_dict\
           else local_dict["functional"]


  def _double_trouble(self):
    """ Returns true, if non-spin polarized or Kammer calculations. """
    from numpy.linalg import norm
    from . import soH
    if self.solo().escan.nbstates  ==   1: return False
    if self.solo().escan.potential != soH: return True
    return norm(self.solo().escan.kpoint) < 1e-12


  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def eigenvalues(self):
    """ Greps eigenvalues from OUTCAR. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    import quantities as pq
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL eigen energies, in eV") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find("*********************************") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    result =  array(result, dtype="float64") if not self._double_trouble()\
              else array([result[i/2] for i in range(2*len(result))], dtype="float64")
    return result * pq.eV

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

    return array(result, dtype="float64") if not self._double_trouble()\
           else array([result[i/2] for i in range(2*len(result))], dtype="float64")

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nnodes(self):
    """ Greps eigenvalue convergence errors from OUTCAR. """
    from os.path import exists, join
    from numpy import array
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
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm), \
            gtor_fourrier( self.raw_gwfns[self.inverse_indices,:,:],\
                           self.rvectors, self.gvectors, self.comm )
        for i, eig in enumerate(self.eigenvalues):
          if i%2 == 0:
            rwfn = rWavefunction( self.comm, i, eig, self._raw_rwfns[0][:,i/2,0],\
                                  self._raw_rwfns[0][:,i/2,1])
          else: 
            rwfn = rWavefunction(self.comm, i, eig, -self._raw_rwfns[1][:,i/2,1].conjugate(),\
                                 self._raw_rwfns[1][:,i/2,0].conjugate())
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
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm), \
            gtor_fourrier( self.raw_gwfns[self.inverse_indices,:,:],\
                           self.rvectors, self.gvectors, self.comm )
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[i%2][:,i/2,0]) )
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i,0]) )
    return result

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
          `lada.physics.reduced_reciprocal_au`.
        - a 3 by x matrix with each row a R-vector in atomic units.
        - one-dimensional array of real coefficients to smooth higher energy G-vectors.
        - one-dimensional array of integer indices to map G-vectors to -G.
    """
    from os import remove
    from os.path import exists, join
    from numpy import sqrt
    from numpy.linalg import norm, det
    from boost.mpi import world
    from quantities import angstrom, pi
    from ..opt import redirect
    from ..opt.changedir import Changedir
    from ..physics import a0, reduced_reciprocal_au
    from ._escan import read_wavefunctions
    from . import soH

    comm = self.comm if self.comm != None else world
    assert self.success
    assert self.nnodes == comm.size, \
           RuntimeError("Must read wavefunctions with as many nodes as they were written to disk.")
    with redirect(fout="") as streams:
      with Changedir(self.directory, comm=self.comm) as directory:
        assert exists(self.escan.WAVECAR),\
               IOError("%s does not exist." % (join(self.directory, self.escan.WAVECAR)))
        self.escan._write_incar(self.comm, self.structure)
        nbstates = self.escan.nbstates if self.escan.potential == soH and norm(self.escan.kpoint)\
                   else self.escan.nbstates / 2
        result = read_wavefunctions(self.escan, range(nbstates), comm)
        remove(self.escan._INCAR + "." + str(world.rank))

    cell = self.structure.cell * self.structure.scale * angstrom
    normalization = det(cell.rescale(a0)) 
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
  def _wavefunction_path(self): return self.solo().escan.WAVECAR

  @property
  def is_spinor(self):
    """ True if wavefunction is a spinor. """
    from . import soH
    return self.escan.potential == soH

  @property
  def is_krammer(self):
    """ True if wavefunction is a spinor. """
    from numpy.linalg import norm
    from . import soH
    return norm(self.escan.kpoint) < 1e-12


  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. `solo` returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    from copy import deepcopy
    
    if self.comm == None: return self
    copy = Extract(self.directory, comm = None)
    copy.OUTCAR = self.OUTCAR
    copy._vffout = self._vffout.solo()
    return copy

  def __repr__(self):
    from os.path import relpath
    return "%s(\"%s\")" % (self.__class__.__name__, relpath(self.directory))

  def __getattr__(self, name):
    """ Passes on public attributes to vff extractor, then to escan functional. """
    if name[0] != '_':
      if hasattr(self._vffout, name): return getattr(self._vffout, name)
      elif self.success: 
        if hasattr(self.escan, name): return getattr(self.escan, name)
    raise AttributeError("Unknown attribute %s" % (name))

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
    if self.success: result.extend( [u for u in dir(self.escan) if u[0] != "_"] )
    return list( set(result) - exclude )

  def __getstate__(self):
    from os.path import relpath
    d = self.__dict__.copy()
    if "comm" in d: del d["comm"]
    if "directory" in d: d["directory"] = relpath(d["directory"])
    return d
  def __setstate__(self, arg):
    self.__dict__.update(arg)
    self.comm = None


