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

""" Submodule to compute bandgaps with escan. """
__docformat__  = 'restructuredtext en'
from ..opt.decorators import make_cached, FileCache
from _extract import Extract as _ExtractE
from .functional import Escan

def extract(outdir=".", comm = None):
  """ Gets extraction object from directory structure. 
  
      Checks first for all-electron calculation in outdir/AE. If it does not
      exist or the calculation is unsucessful, then checks for VBM, CBM
      directories, and checks those calculations. If unsucessfull, returns a
      fake extraction object with success set to False. Otherwise, returns
      successfull extraction object.
  """
  from os.path import exists, join
  from ..mpi import Communicator
  from ..opt import RelativeDirectory

  comm = Communicator(comm)
  outdir = RelativeDirectory(outdir).path

  paths = join(outdir, "AE"), join(outdir, "VBM"), join(outdir, "CBM")
  exists_paths = comm.broadcast([exists(p) for p in paths] if comm.is_root else None)
  if exists_paths[0]: 
    result = ExtractAE( _ExtractE(paths[0], comm = comm) )
    return result
  elif exists_paths[1] and exists_paths[2]:
    result = ExtractRefs( _ExtractE(paths[1], comm = comm),\
                          _ExtractE(paths[2], comm = comm),
                          _ExtractE(outdir, comm = comm) )
    return result
  class NoBandGap(object): 
    @property
    def success(self): return False
    @property
    def directory(self): return outdir
    def __str__(self): return "NoBandGap({0})".format(outdir)
    def __repr__(self): return "NoBandGap({0})".format(outdir)
  return NoBandGap()

def bandgap(escan, structure, outdir=None, references=None, n=5, overwrite = False, **kwargs):
  """ Computes bandgap of a structure with a given escan functional. 
  
      The band-gap is computed using an all-electron method (if references=None
      in argument), or a folded spectrum method. The latter expects two
      references: one slightly above the vbm, and the other slightly below the
      cbm. It performs two folded spectrum calculations, one for each method,
      using the functional and keyword arguments given on input. If the
      references are not placed accurately (eg slightly above the CBM and
      slightly below the VBM), the algorithm applies some heuristics to try an
      determine a better set of references. The calculations then repeat for a
      maximum of n times. Beyond that, or if the references vs eigenvalues
      cannot be made sense of, an electron-calculation is performed.
  """
  from os import getcwd
  from os.path import abspath
  from copy import deepcopy
  from ..mpi import Communicator

  assert "do_escan" not in kwargs,\
         ValueError("\"do_escan\" is not an admissible argument of bandgap.")

  escan = deepcopy(escan)
         
  if outdir is None: outdir = getcwd()
  outdir    = abspath(outdir)
  overlap_factor = kwargs.pop("overlap_factor", 10e0)

  comm = Communicator(kwargs.pop("comm", None), with_world=True)
  if not overwrite:  # check for previous results.
    result = extract(outdir, comm)
    if result.success: return result

  
  kwargs["overwrite"] = overwrite
  kwargs["comm"] = comm
  return _band_gap_ae_impl(escan, structure, outdir, **kwargs) if references is None\
         else _band_gap_refs_impl(escan, structure, outdir, references, n, \
                                  overlap_factor=overlap_factor, **kwargs) 

class ExtractAE(object):
  """ Band-gap extraction class. """
  is_ae = True
  """ This was an all-electron bandgap calculation. """
  def __init__(self, extract):
    """ Initializes extraction object for all-electron bandgaps. """
    self.__dict__['_extract'] = extract
    super(ExtractAE, self).__init__()

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  def vbm(self):
    """ Greps VBM from calculations. """
    return self._vbm_cbm[0]

  @property
  def cbm(self):
    """ Greps CBM from calculations. """
    return self._vbm_cbm[1]

  @property
  @make_cached
  def _vbm_cbm(self):
    """ Gets vbm and cbm. """
    from quantities import eV
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    n = len(self.structure.atoms) * 4
    a, b = eigenvalues[n-1], eigenvalues[n] 
    if hasattr(a, "units"): a.units = eV
    else: a = a * eV
    if hasattr(b, "units"): b.units = eV
    else: b = b * eV
    return a, b


  def oscillator_strength(self, degeneracy=1e-3, attenuate=False):
    """ Computes oscillator strength between vbm and cbm. """
    from numpy import abs, dot
    from ..physics import electronic_mass, h_bar

    units = 2e0/3e0 * h_bar**2 / electronic_mass
    result, nstates = None, 0
    for eigA, eigB, dipole in self.dipole(degeneracy, attenuate):
      if abs(eigA - self.cbm) > degeneracy: continue
      if abs(eigB - self.vbm) > degeneracy: continue
      nstates += 1
      dme = dot(dipole, dipole.conjugate()).real * dipole.units * dipole.units
      if result is None: result = dme / (eigA - eigB) 
      else: result += dme / (eigA - eigB) 
    return (units * result).simplified, nstates


  def dipole(self, degeneracy=-1e0, attenuate=False):
    """ Computes dipole matrix element between vbm and cbm. """
    # gets result, possibly from cache file.
    d2, a2, result = self._dipole(degeneracy, attenuate)

    uncache  = degeneracy < 0e0 and d2 >= 0e0
    if d2 > 0e0 and degeneracy > 0e0: 
      uncache |= abs(d2 - degeneracy) >= min(d2, degeneracy) and d2 < degeneracy
    uncache |= a2 != attenuate
    if uncache: 
      from os.path import join
      from os import remove
      remove(join(self._extract.directory, "DIPOLECAR"))
      return self._dipole(degeneracy, attenuate)[-1]
    return result
    

  @FileCache('DIPOLECAR')
  def _dipole(self, degeneracy=-1e0, attenuate=False):
    """ Computes dipole matrix element between vbm and cbm. 
    
        This routine caches results in a file. The routine above should check
        that the arguments are the same.
    """
    from ..physics import a0
    result, gvectors = [], self.gvectors.rescale(1./a0)
    for wfnA in self.gwfns:
      if degeneracy >= 0e0 and abs(wfnA.eigenvalue - self.cbm) > degeneracy: continue
      for wfnB in self.gwfns:
        if degeneracy >= 0e0 and abs(wfnB.eigenvalue - self.vbm) > degeneracy: continue
        result.append( (wfnA.eigenvalue, wfnB.eigenvalue,
                        wfnA.braket(gvectors, wfnB, attenuate=attenuate)) )
    return degeneracy, attenuate, result

  def __copy__(self):
    """ Returns a shallow copy of this object. """
    result = self.__class__(self)
    return result

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        :kwarg dipolecar: Include dipole moment file.
        :type dipolecar: bool

        Parameters are passed on to internal escan calculation.
    """
    if kwargs.get('dipolecar', False): 
      from os.path import exists, join
      if exists(join(self._extract.directory, 'DIPOLECAR')):
        yield join(self._extract.directory, 'DIPOLECAR')
    for file in self._extract.iterfiles(**kwargs): yield file

  @property
  def directory(self):
    """ Directory containing the results. """
    from os.path import dirname, realpath
    return dirname(realpath(self._extract.directory))


  def __getattr__(self, name):
    """ Passes on result to extract object. """
    return getattr(self._extract, name)

  def __dir__(self):
    """ All public attributes in this object. """
    return   [u for u in dir(self.__class__) if u[0] != '_'] \
           + [u for u in dir(self._extract) if u[0] != '_']
  
def _band_gap_ae_impl(escan, structure, outdir, **kwargs):
  """ Computes bandgap of a structure using all-electron method. """
  from os.path import join
  
  if "eref" in kwargs:
    assert kwargs["eref"] is None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  outdir = join(outdir, "AE")
  nbstates = kwargs.pop("nbstates", escan.nbstates)
  if nbstates is None: nbstates = 4
  if nbstates == 0: nbstates = 4
  if hasattr(nbstates, '__getitem__'):
    assert len(nbstates) > 1, ValueError("Not sure what nbstates is.")
    nbstates = max(nbstates)
  nbstates = nbstates  + len(structure.atoms) * 4
  extract = escan( structure, outdir = outdir, eref = None,\
                   nbstates = nbstates, **kwargs)
  return ExtractAE(extract)

class ExtractRefs(object):
  """ Band-gap extraction class for folded spectrum 2-ref method. """
  is_ae = False
  """ This was not an all-electron bandgap calculation. """
  def __init__(self, vbm_out, cbm_out, vff_out):
    super(ExtractRefs, self).__init__()
    self.extract_vbm = vbm_out
    """ VBM extraction method. """
    self.extract_cbm = cbm_out
    """ CBM extraction method. """
    self.extract_vff = vff_out
    """ VFF extraction method. """
    self.comm = self.comm # makes sure all share the same communicator.

  @property
  def comm(self):
    """ Communicator over which to sync output. """
    return self.extract_vbm.comm
  @comm.setter
  def comm(self, value):
    self.extract_vbm.comm = value
    self.extract_cbm.comm = value
    self.extract_vff.comm = value

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  def _raw(self):
    from numpy import array
    from quantities import eV
    vbm_eigs = self.extract_vbm.eigenvalues.copy()
    cbm_eigs = self.extract_cbm.eigenvalues.copy()

    vbm_ref = self.extract_vbm.functional.eref
    cbm_ref = self.extract_cbm.functional.eref
    above_refs =       [u for u in cbm_eigs if u >= cbm_ref]
    above_refs.extend( [u for u in vbm_eigs if u >= cbm_ref])
    above_refs = array([u for u in set(above_refs)] );
    above_refs.sort()
    below_refs =       [u for u in cbm_eigs if u <= vbm_ref]
    below_refs.extend( [u for u in vbm_eigs if u <= vbm_ref])
    below_refs = array([u for u in set(below_refs)] );
    below_refs.sort()

    a, b = below_refs[-1], above_refs[0]
    if hasattr(a, "units"): a.units = eV
    else: a = a * eV
    if hasattr(b, "units"): b.units = eV
    else: b = b * eV
    return a, b

  @property
  def vbm(self):
    """ Greps VBM from calculations. """
    return self._raw[0]

  @property
  def cbm(self):
    """ Greps CBM from calculations. """
    return self._raw[1]

  def oscillator_strength(self, degeneracy=1e-3, attenuate=False):
    """ Computes oscillator strength between vbm and cbm. """
    from numpy import abs, dot
    from ..physics import electronic_mass, h_bar

    units = 2e0/3e0 * h_bar**2 / electronic_mass
    result, nstates = None, 0
    for eigA, eigB, dipole in self.dipole(degeneracy, attenuate=attenuate):
      if abs(eigA - self.cbm) > degeneracy: continue
      if abs(eigB - self.vbm) > degeneracy: continue
      nstates += 1
      dme = dot(dipole, dipole.conjugate()).real * dipole.units * dipole.units
      if result is None: result = dme / (eigA - eigB) 
      else: result += dme / (eigA - eigB) 
    return (units * result).simplified, nstates
  
  def dipole(self, degeneracy=-1e0, attenuate=False):
    """ Computes dipole matrix element between vbm and cbm. """
    # gets result, possibly from cache file.
    try:  d2, a2, result = self._dipole(degeneracy, attenuate)
    except ValueError: # there might be an issue with older DIPOLECAR.
      result = self._dipole(degeneracy, attenuate)
      uncache = False
    else:
      uncache  = degeneracy < 0e0 and d2 >= 0e0
      if d2 > 0e0 and degeneracy > 0e0: 
        uncache |= abs(d2 - degeneracy) >= min(d2, degeneracy) and d2 < degeneracy
      uncache |= a2 != attenuate
    if uncache: 
      from os.path import join
      from os import remove
      remove(join(self.directory, "DIPOLECAR"))
      return self._dipole(degeneracy, attenuate)[-1]
    return result
    
  @FileCache('DIPOLECAR')
  def _dipole(self, degeneracy=-1e0, attenuate=False):
    """ Computes dipole matrix element between vbm and cbm. """
    from numpy import all, abs
    from ..physics import a0

    assert self.extract_vbm.gvectors.shape == self.extract_cbm.gvectors.shape
    assert all( abs(self.extract_vbm.gvectors - self.extract_cbm.gvectors) < 1e-12 )
    result, gvectors = [], self.extract_cbm.gvectors.rescale(1./a0)
    for wfnA in self.extract_cbm.gwfns:
      if degeneracy >= 0e0 and abs(wfnA.eigenvalue - self.cbm) > degeneracy: continue
      for wfnB in self.extract_vbm.gwfns:
        if degeneracy >= 0e0 and abs(wfnB.eigenvalue - self.vbm) > degeneracy: continue
        result.append( (wfnA.eigenvalue, wfnB.eigenvalue, 
                       wfnA.braket(gvectors, wfnB, attenuate=attenuate)) )
    return degeneracy, attenuate, result

  @property
  def success(self):
    """ True if all calculations were successfull. """
    return self.extract_vff.success and self.extract_vbm.success and self.extract_cbm.success

  def __getattr__(self, name): 
    """ Sends to vff output. """
    if hasattr(self.extract_vff, name): return getattr(self.extract_vff, name)
    raise AttributeError("Unknown attribute %s." % (name))

  def __dir__(self):
    """ Returns list of attributes. 

        Since __getattr__ is modified to include vff extraction data, __dir__
        should be modified as well. This makes for better ipython
        (tab-completion) integration. 
    """
    result = [u for u in dir(self.__class__) if u[0] != '_'] 
    result.extend([u for u in self.__dict__.keys() if u[0] != '_'])
    result.extend(dir(self.extract_vff))
    return result

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        :kwarg dipolecar: Include dipole moment file.
        :type dipolecar: bool

        Parameters are passed on to internal escan calculation.
    """
    if kwargs.get('dipolecar', False): 
      from os.path import exists, join
      if exists(join(self.directory, 'DIPOLECAR')): yield join(self.directory, 'DIPOLECAR')
    for file in self.extract_vbm.iterfiles(**kwargs): yield file
    for file in self.extract_cbm.iterfiles(**kwargs): yield file
    try: extract = self.extract_vbm.__class__( dirname(self.directory), 
                                               comm  = self.comm,
                                               escan = self.extract_vbm.functional )
    except: pass
    else:
      for file in extract.iterfiles(**kwargs): yield file

def _band_gap_refs_impl( escan, structure, outdir, references, n=5,\
                         overlap_factor=10e0, **kwargs):
  """ Computes band-gap using two references. """
  from os.path import join
  from numpy import array, argmax
  
  # check/correct input arguments
  if "overwrite" in kwargs: del kwargs["overwrite"]
  nbstates = kwargs.pop("nbstates", escan.nbstates)
  if nbstates < 2: nbstates = 2
  if "eref" in kwargs:
    assert kwargs["eref"] is None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  vffrun = kwargs.pop("vffrun", escan.vffrun)
  genpotrun = kwargs.pop("genpotrun", escan.genpotrun)
  
  assert len(references) == 2, ValueError("Expected 2-tuple for argument \"references\".")
  vbm_ref, cbm_ref = references
  if vbm_ref > cbm_ref: cbm_ref, vbm_ref = references

  # first computes vff and genpot unless given.
  if genpotrun is None or vffrun is None: 
    vffout = escan( structure, outdir=outdir, do_escan=False, genpotrun=genpotrun,\
                    vffrun=vffrun, **kwargs )
    if genpotrun is None: genpotrun = vffout
    if vffrun is None: vffrun = vffout
  else: vffout = vffrun

  iter, continue_loop = 0, True
  recompute = [True, True]
  while iter < n and continue_loop:
    # computes vbm
    if recompute[0]:
      try: n = nbstates[0]
      except: n = nbstates
      vbm_out = escan\
                (
                  structure, outdir=join(outdir,"VBM"), \
                  eref=vbm_ref, overwrite=True, vffrun=vffrun,\
                  genpotrun=genpotrun, nbstates=n, **kwargs 
                )
      vbm_eigs = vbm_out.eigenvalues.copy()
    # computes cbm
    if recompute[1]:
      try: n = nbstates[1]
      except: n = nbstates
      cbm_out = escan\
                (
                  structure, outdir=join(outdir, "CBM"), \
                  eref=cbm_ref, overwrite=True, vffrun=vffrun,
                  genpotrun=genpotrun, nbstates=n, **kwargs
                )
      cbm_eigs = cbm_out.eigenvalues.copy()
    recompute = [False, False] # by default, does not recompute
  
    below_refs =       [u for u in cbm_eigs if u < vbm_ref]
    below_refs.extend( [u for u in vbm_eigs if u < vbm_ref])
    below_refs = array([u for u in set(below_refs)] )
    below_refs.sort()
    between_refs =       [u for u in cbm_eigs if u >= vbm_ref and u < cbm_ref]
    between_refs.extend( [u for u in vbm_eigs if u >= vbm_ref and u < cbm_ref])
    between_refs = array([u for u in set(between_refs)] );
    between_refs.sort()
    above_refs =       [u for u in cbm_eigs if u >= cbm_ref]
    above_refs.extend( [u for u in vbm_eigs if u >= cbm_ref])
    above_refs = array([u for u in set(above_refs)] )
    above_refs.sort()
  
    if between_refs.size == 0: # no eigenvalues between the references.
      if below_refs.size > 0 and above_refs.size > 0: break # sole case where break is allowed.
      continue_loop = False; continue # got to all electron calculation.

    # there are eigenvalues between the references. Determines the largest "gap"
    a = [ vbm_ref ]; a.extend(u for u in between_refs.flat); a.append(cbm_ref)
    gap_index = argmax(array(a[1:]) - array(a[:-1])) # computes all differences.
    deltas = [a[gap_index] - a[0], a[gap_index+1] - a[gap_index], a[-1] - a[gap_index+1]]
    # Check pathological case where vbm and cbm give essentially same eigenvalues.
    if gap_index == 0 and below_refs.size == 0:
      if deltas[1] > overlap_factor * deltas[2]: vbm_ref -= deltas[1] * 0.95; recompute[0] = False
      else: continue_loop = False;
      continue 
    if gap_index == len(a)-1 and above_refs.size == 0:
      if deltas[1] > overlap_factor * deltas[1]: cbm_ref += deltas[1] * 0.95; recompute[1] = False
      else: continue_loop = False;
      continue #
    # Computes actual gap.
    if gap_index == 0: deltas[1] += vbm_ref - below_refs[-1]
    if gap_index == len(a)-1: deltas[1] += above_refs[0] - cbm_ref 
    # case where no gap can be truly determined.
    if not (deltas[1] > overlap_factor * deltas[0] and deltas[1] > overlap_factor * deltas[2]):
      continue_loop = False; continue # go to all electron calculation.

    # finally, recomputation case. Sets reference to best possible values.
    vbm_ref = below_refs[-1] if gap_index == 0        else a[gap_index]
    cbm_ref = above_refs[0]  if gap_index == len(a)-1 else a[gap_index+1]
    vbm_ref +=  deltas[1] * 0.3 
    cbm_ref -=  deltas[1] * 0.3 
    recompute = [True, True]
    iter += 1
  else: # ran through all iterations and failed.
    return _band_gap_ae_impl(escan, structure, outdir, **kwargs)
  return ExtractRefs(vbm_out, cbm_out, vffout)

class Functional(Escan): 
  """ Bandgap functional.
  
      Computes bandgap using either full diagonalization or two folded spectrum
      calculations. Performs some checking so that a non-zero band-gap is
      discored. see `bandgap` method.
  """
  Extract = staticmethod(extract)
  def __init__(self, *args, **kwargs):
    """ Initializes an LDOS functional. 
    
        :param args: Any argument that works for `Escan`.
        :param kwargs: Any keyword argument that works for `Escan`.
    """
    self.references = kwargs.pop('references', None)
    """ References for folded spectrum calculation, or None for full diagonalization. 

        In the first case, this is a tuple giving the two reference energies
        (CBM and VBM) in eV.
    """ 
    self.n = kwargs.pop('n', 5)
    """ Maximum number of trial folded-spectrum calculations. 

        Beyond this, resorts to full diagonalization.
    """ 
    self.dipole = kwargs.pop('dipole', False)
    """ Whether or not to compute dipole elements. """
    escan_copy = kwargs.pop('escan', None)
    super(Functional, self).__init__(**kwargs)

    # copies parent functional.
    if escan_copy is not None:
      from copy import deepcopy
      self.__dict__.update(deepcopy(escan_copy.__dict__))

  def __call__(self, structure, outdir=None, comm=None, **kwargs):
    """ Computes band-gap. 

        Parameters are passed on to `bandgap` method.
    """
    if '_computing' in self.__dict__:
      return super(Functional, self).__call__(structure, outdir, comm, **kwargs)

    self._computing = True
    try: 
      if 'references' not in kwargs: kwargs['references'] = self.references
      if 'n' not in kwargs: kwargs['n'] = self.n
      dipole = kwargs.pop('dipole', self.dipole)
      escan = Escan()
      escan.__dict__.update(self.__dict__)
      result = bandgap(escan, structure, outdir, comm=comm, **kwargs)
      if dipole: result.dipole()
      return result
    finally: del self._computing
