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

""" Functional to compute effective mass. """
__docformat__ = "restructuredtext en"
__all__ = ["Functional", "Extract"]
from ..opt.decorators import make_cached
from .kescan import Extract as KExtract, KEscan

class Extract(KExtract):
  """ Extraction object for effective mass functional. """
  def __init__(self, directory=None, comm=None, bandgap=None, **kwargs):
    """ Initializes the extraction object. """
    super(Extract, self).__init__(directory, comm, **kwargs)
    if bandgap is not None: self._extract_bg = bandgap
  
  @property
  def type(self): 
    """ Whether electronic or hole effective masses are computed. """
    return self.functional.type

  @property
  def direction(self): 
    """ Direction for which effective mass is computed. """
    return self.functional.direction

  @property
  def center(self): 
    """ K-point at which effective mass is computed. """
    return self.functional.center 

  @property
  def extract_bg(self, *args, **kwargs):
    """ Band-gap extractor object. """
    if hasattr(self, "_extract_bg"): return self._extract_bg
    from . import extract_bg as ebg
    e = ebg(self.directory, self.comm)
    assert e.success, RuntimeError("Could not extract bandgap from calculations at {0}.".format(self.directory))
    return e

  @property 
  def cbm(self):
    """ Conduction Band Minimum at kpoint for which effective mass is computed. """
    return self.extract_bg.cbm 

  @property 
  def vbm(self):
    """ Valence Band Minimum at kpoint for which effective mass is computed. """
    return self.extract_bg.vbm 

  @property 
  def bandgap(self):
    """ Band-gap at kpoint for which effective mass is computed. """
    return self.extract_bg.bandgap
  
  def dipole(self, *args, **kwargs):
    """ Computes dipole element between vbm and cbm. """
    return self.extract_bg.dipole(*args, **kwargs)

  def oscillator_strength(self, *args, **kwargs):
    """ Computes dipole element between vbm and cbm. """
    return self.extract_bg.dipole(*args, **kwargs)

  @property 
  def _is_array(self):
    """ True if direction is array of directions. """
    from numpy import array
    direction = array(self.functional.direction)
    return direction.ndim == 2
  @property
  def _nbpoints(self):
    """ Meaningful number of points. """
    return max(3, self.functional.nbpoints)

  @property
  @make_cached
  def derivatives(self): 
    """ Effective mass in given direction. 
    
        Returned in units of the electron's mass at rest. Sign difference
        between electron and holes is accounted for.
    """
    from numpy import sort, array
    from numpy.linalg import lstsq as np_lstsq
    from quantities import hartree
    lstsq = self.functional.lstsq
    if lstsq is None: lstsq = np_lstsq

    # constructs matrix of eigenvalues with cbm or vbm only.
    eigs = []
    eigenvalues = KExtract.eigenvalues.__get__(self).rescale(hartree)
    if self.type == "e": 
      for e in eigenvalues:
        eigs.append( e[e > (self.vbm + 0.5 * self.bandgap)] )
    elif self.type == "h":
      for e in eigenvalues:
        eigs.append( e[e < (self.vbm + 0.5 * self.bandgap)] )
    else: raise RuntimeError("type ({0}) is neither \"e\" nor \"h\".")
    m = min([len(e) for e in eigs])
    eigs = array([ e[:m] for e in eigs ]) * hartree
    measurements = sort(eigs, axis=1) 
    # makes sure that parameters are constructed.
    kpoints = self.functional.kpoints
    for k in kpoints(self.input_structure, self.structure): continue
    # depending on format of input directions, returns derivatives.
    if self._is_array:
      result = []
      parameters = kpoints.parameters
      for i in xrange(kpoints.directions.shape[0]):
        j, k = i * self._nbpoints, (i+1) * self._nbpoints
        result.append(lstsq(parameters[i,:,:], measurements[j:k, :])[0])
      return array(result)
    return lstsq(kpoints.parameters, measurements)[0]

  @property
  def eigenvalues(self): 
    """ Eigenvalues at k-point where effective masses are computed. 
    
        These eigenvalues are the result of the least-square-fit.
    """
    from quantities import eV, hartree
    if self._is_array:
      result = (self.derivatives[:,0,:] * hartree).rescale(eV)
      return result if self.type == "e" else result[:,::-1]
    else:
      result = (self.derivatives[0,:] * hartree).rescale(eV) 
      return result if self.type == "e" else result[::-1]

  @property
  def mass(self): 
    """ Effective mass in given direction. 
    
        Returned in units of the electron's mass at rest. Sign difference
        between electron and holes is accounted for.
    """
    assert self.type == "e" or self.type == "h",\
           RuntimeError("Unknown type {0}.".format(self.type))
    if self._is_array:
      return 1e0/self.derivatives[:,2,:] if self.type == "e" else -1e0/self.derivatives[:,2,::-1]
    else:
      return 1e0/self.derivatives[2,:] if self.type == "e" else -1e0/self.derivatives[2,::-1]
  
  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        Parameters are passed on to internal escan calculations.
    """
    from itertools import chain
    return chain( self.extract_bg.iterfiles(**kwargs), 
                  super(Extract, self).iterfiles(**kwargs) )

class Functional(KEscan):
  """ Effective mass functional. 
  
      The effective mass is computed for a given set of bands, for the VBM, or
      the CBM, in a  specific direction. 
  """
  Extract = Extract
  """ Extraction object for the effective-mass functional. """
  def __init__( self, direction=None, nbpoints=0, stepsize=1e-2, \
                center=None, lstsq=None, **kwargs ):
    """ Computes effective mass for a given direction.
    
        :Parameters:
          type : "e" or "h"
            Whether to compute electronic or effective mass.
          direction : 3-tuple 
            direction for which to compute effective mass.
            Defaults to (0,0,1).
          nbpoints : int
            Number of points with wich to compute derivatives.
            Should be at least order + 1. Default = order + 1. 
          stepsize : float
            Distance between interpolation points. Default = 1e-2.
            Units of ``2|pi|/a``, with ``a=structure.scale``.
          center : 3-tuple
            k-point where to take derivative. Units of ``2|pi|/a``, with
            ``a=structure.scale``.
            Defaults to (0,0,0).
          lstsq 
            Linear least square method. The first two parameters should
            be same as numpy.linalg.lstsq. Other parameters can be passed as extra
            parameters. Defaults to numpy.linalg.lstsq.
          kwargs 
            Extra parameters which are passed on first to escan (if escan
            object has an attribute with the same name), then to the linear least
            square fit method. Note that this will *not* change the external escan
            object.  This function is stateless. 
    
        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    from numpy import array
    if center is None: center = (0, 0, 0)
    if direction is None: direction = (0, 0, 1)

    self.type      = "e"
    """ Whethe to compute electronic or hole effective masses. """
    self.direction = array(direction, dtype="float64")
    """ Direction for which to compute effective mass. """
    self.nbpoints = nbpoints
    """ Number of points with which to perform least-square fit. Defaults to 3. """
    self.stepsize = stepsize
    """ Distance between interpolation points. Default = 1e-2.

        Units of ``2|pi|/a``, with ``a=structure.scale``.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.lstsq = lstsq
    """ Least square fit method to use when computing effective mass. 
    
        If None, will default ot numpy.linalg.lstsq. Otherwise, it should be a
        pickleable-callable, or suffer the consequences.
    """
    super(Functional, self).__init__(**kwargs)
    
    # gotta be after super call, otherwise overwritten.
    self.center = center
    """ k-point for which to compute effective mass. """
    self.__dict__.pop('kpoints', None)

  @property
  def kpoints(self):
    """ KPoint derived class for effective mass. """
    from numpy import array
    from .derivatives import ReducedDDPoints, ReducedChainedDDPoints

    kwargs = { "direction": array(self.direction, dtype="float64"),
               "center": self.center,
               "order": 2,
               "nbpoints": max(self.nbpoints, 3),
               "stepsize": self.stepsize,
               "relax": self.do_relax_kpoint }

    return ReducedChainedDDPoints(**kwargs) if kwargs['direction'].ndim == 2 \
           else ReducedDDPoints(**kwargs)
  @kpoints.setter
  def kpoints(self, value): pass
    
  @property 
  def center(self):
    """ k-point for which to compute effective mass. 
    
        This object is owned by the functional, i.e. a copy of ``something`` is
        created when setting ``self.direction = something``.
    """
    return self.kpoint
  @center.setter
  def center(self, value):
    from numpy import array
    self.kpoint = array(value, dtype="float64")

  def __call__(self, structure, outdir=None, comm=None, bandgap=None, references=None, n=5, **kwargs):
    """ Computes effective mass.

        :Parameters: 
          structure : `crystal.Structure`
            Structure for which to compute effective mass.
          outdir : str or None
            Directory where to save results. If None, will create results are
            stored in current directory.
          comm : None or `mpi.Communicator`
            MPI Communicator holding processes with which to perform
            calculations.
          bandgap 
            If given, should be the return from a bandgap calculation.
          references
            For bandgap calculation, as necessary. See BandGap.
          n
            For bandgap calculation, as necessary. See BandGap.
        Parameters are passed on to `dervatives.reciprocal` method.
    """
    from ..mpi import Communicator
    from ._bandgap import Functional as Bandgap
    if '_computing' in self.__dict__:
      return super(Functional, self).__call__(structure, outdir, comm, **kwargs)

    if comm is None: comm = Communicator(comm, with_world=True)

    # copy functional with current type. noadd keyword makes sure that only
    # known attributes are added.
    this = self.copy(noadd=True, **kwargs)

    # symlinks files from bandgap calculations.
    if bandgap is not None: this._link_bg_files(bandgap, outdir, comm)

    # computes bandgap.
    bandgap_func = Bandgap(escan=this, references=references, n=n)
    bandgap = bandgap_func(structure, outdir, comm, **kwargs)
    assert bandgap.success, RuntimeError("Could not compute bandgap.")
    
    # then figures out appropriate reference.
    if this.type == "e":   eref = bandgap.cbm - bandgap.bandgap / 4e0
    elif this.type == "h": eref = bandgap.vbm + bandgap.bandgap / 4e0
    else: raise ValueError("Unknown mass type {0}.".format(this.type))
    
    # performs calculation.
    kout = super(Functional, this).__call__(structure, outdir, comm=comm, eref=eref, **kwargs)
    # Saves FUNCCAR.
    if comm.is_root:
      from os.path import join
      from cPickle import dump
      with open(join(kout.directory, this._FUNCCAR), "w") as file: dump(this, file)

    # Effective mass extractor.
    return self.Extract(outdir, comm, unreduce=True, bandgap=bandgap)


  def _link_bg_files(self, bandgap, outdir, comm):
    """ Creates link to band-gap files in output directory. """
    from os import symlink
    from os.path import join, exists, relpath, dirname, basename
    from ..opt import Changedir

    if bandgap is None: return
    assert bandgap.success,\
           RuntimeError( "Input bandgap calculations at {0} were not successfull."\
                         .format(bandgap.directory) )
    result = self.Extract(outdir, comm, unreduce=True, bandgap=bandgap)
    directory = bandgap.directory
    if bandgap.is_ae: directory = dirname(directory)
    if comm.is_root: 
      for file in bandgap.iterfiles(): 
        name = join(result.directory, relpath(file, bandgap.directory))
        if not exists(name): # tries and links to earliest directory.
          alldirs = relpath(file, bandgap.directory).split('/')
          src, dest = bandgap.directory, result.directory
          for now in alldirs:
            src, dest = join(src, now), join(dest, now)
            with Changedir(dirname(dest)) as cwd:
              if not exists(basename(dest)):
                symlink(relpath(src, dirname(dest)), basename(dest))
                break
