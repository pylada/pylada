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

""" Defines Local Density of States.

    Both a functional (`Functional`) and a function (`ldos`) are defined which
    compute the Local Density of State.  The LDOS is obtaine as a Fourier
    transform of the wavefunctions to the input real-space positions. In this
    way, positions are not limited to the FFT mesh.

    The functional returns an extraction object which caches the LDOS to file.
"""
__docformat__ = "restructuredtext en"
__all__ = ['ldos', 'Extract', 'Functional']

from .kescan import KEscan, Extract as KExtract
from ..opt import make_cached, FileCache

class _ldosfunc(object):
  """ Local density of states for a given set of positions within a given structure. """
  def __init__(self, eigenvalues, rs, volume):
    """ Initializes a local density of state functor. 
    
        :Parameters:
          eigenvalues : numpy array
            Vector of eigenvalues. Can be signed with a unit, or without, in
            which case eV is assumed.
          rs : numpy array
            Matrix of densities per real-space position(row) and per band(column). 
          volume : scalar
            Volume of the brillouin zone.
    """
    from numpy import sqrt, pi
    from quantities import eV
    self.eigenvalues = eigenvalues.copy()
    """ Vector of eigenvalues. """
    self.rs = rs.copy()
    """ Matrix of densities per real-space position(row) and per band(column). """
    if not hasattr(self.eigenvalues, 'rescale'): self.eigenvalues *= eV 
    else: self.eigenvalues = self.eigenvalues.rescale(eV)
    self.normalization = float(volume) / float(self.eigenvalues.shape[0]) / sqrt(pi)
    """ Normalization for LDOS.
    
        brillouin zone volume / number of k-points / sqrt(|pi|). 

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """

  def __call__(self, energy, sigma=0.1):
    """ Calls smearing function over densities.
    
        :Parameters: 
          energy : float or scalar array
            Energy at which to compute local density of states. This can be real
            number, in which case it should be in eV, or a numpy scalar with a
            unit (from quantity).
          sigma : float or scalar array
            Width at half maximum of the gaussian smearing. This can be real
            number, in which case it should be in eV, or a numpy scalar with a
            unit (from quantity).
    """
    from numpy import dot, exp, array
    from quantities import eV
    if not hasattr(sigma, 'rescale'): sigma *= eV
    else: sigma = sigma.rescale(eV)
    if not hasattr(energy, 'rescale'): energy = array(energy) * eV
    else: energy = energy.rescale(eV)
    y = array([ [(E - e)/sigma for e in energy] for E in self.eigenvalues])
    return dot(self.rs, exp(-y*y)) * (float(self.normalization) / sigma)


def ldos(extractor, positions, raw=False):
  """ Local density of states at given positions.
  
      :Parameters:
        extractor 
          Output from a `KEscan` or a `Functional`  calculation.
        positions : nx3 array
          Array with positions in real space where to compute LDOS.
        raw : boolean
          Whether to return the raw data or the LDOS itself, i.e. a function of
          the energy.
  """
  from numpy import tensordot, multiply, conjugate, exp, concatenate,\
                    array, rollaxis, sum, add, zeros
  from numpy.linalg import det, inv
  from quantities import angstrom

  assert isinstance(extractor, KExtract),\
         ValueError('extractor argument should be KExtract isntance.')

  extractor = extractor.copy(unreduce=False)
  istr, ostr = extractor.input_structure, extractor.structure
  normalization = 0e0
  perpoint = []
  for n, equivs in enumerate(extractor.functional.kpoints.iter_equivalents(istr, ostr)):

    extract = extractor[n]

    # checks that this proc can return wavefunctions. 
    is_null = extract.raw_gwfns is None
    Neigs, Npos = extract.eigenvalues.shape[0], len(positions)
    extract = extract.copy(comm=extract.comm.split(0 if is_null else 1))
    if is_null: perpoint.append(zeros(shape=(Npos, Neigs)))
    else: 
      # computes all positions including symmetry equivalents.
      # Since we expect fewer real space points than fourrier space points,
      # symmetric equivalents are added to real-space positions. 
      # See for instance "Electronic Structure", Richard M. Martin, first
      # edition, chapter 4 section 5.
      equivs = [u for u in equivs]
      operators = [op.inverse for index, m, k, op in equivs]
      all_positions = array([op(u) * getattr(u, "units", angstrom) for op in operators for u in positions])
      for u in all_positions: 
        if hasattr(u, "units"): u.rescale(angstrom)
      all_positions = array(all_positions) * angstrom
      multiplicities = [m for index, m, k, op in equivs]
      normalization += sum(multiplicities)
  
      # creates array which may include krammer degenerate.
      if extract.is_krammer:
        inverse = conjugate(extract.raw_gwfns[extract.inverse_indices,:,:])
        gwfns = concatenate((extract.raw_gwfns, inverse), axis=1)
      else: gwfns = extract.raw_gwfns
      # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
      # units are not conserved by tensordot, so must do it by hand.
      units = (all_positions.units * extract.gvectors.units).simplified
      v = exp(-1j * tensordot(all_positions, extract.gvectors, ((1),(1))) * units)
      # computes fourrier transform for all wavefunctions simultaneously.
      rspace = tensordot(v, gwfns, ((1),(0)))
      rspace = multiply(rspace, conjugate(rspace)).real
      # Sum over spin channels if necessary.
      if not extract.is_spinor: rspace = rspace[:,:,0]
      else: rspace = rspace[:,:,0] + rspace[:,:,1]
      # Sum degenerate states if necessary.
      if extract.is_krammer:
        assert rspace.shape[1] % 2 == 0
        # reorder array same as eigenvalues.
        result = array([ (rspace[:,i//2] if i % 2 == 0 else rspace[:,i//2+rspace.shape[1]//2])\
                         for i in xrange(rspace.shape[1]) ])
      
      # sum over equivalent kpoints. 
      if abs(multiplicities[0] - 1e0) > 1e-12: rspace[:Npos, :] *= multiplicities[0]
      for j, m in enumerate(multiplicities[1:]):
        if abs(m - 1e0) > 1e-12: rspace[:Npos, :] += m * rspace[(j+1)*Npos:(j+2)*Npos, :]
        else: rspace[:Npos, :] += rspace[(j+1)*Npos:(j+2)*Npos, :]
  
      # append to reduced kpoint ldos list.
      perpoint.append(rspace[:Npos,:].copy())

  # normalize results and concatenate.
  result = rollaxis(array(perpoint), 0,-1) / float(normalization)
  result = extractor.comm.all_reduce(array([u.flatten() for u in result]), add)
  
  return result if raw \
         else _ldosfunc(extractor.eigenvalues, result, det(inv(extractor.structure.cell)))



class Extract(KExtract):
  """ Extraction routine for LDOS. """
  def __init__(self, *args, **kwargs): 
    """ Creates Extraction object. 


        All parameters are passed on to KExtract.__init__, unless
        ``parent`` is present. In that case, ``parent`` should be a KExtract
        object which will be copied. This way, we can add ldos specific
        properties to a KExtract object.
    """
    parent = kwargs.pop('parent', None)
    if parent is not None:
      assert len(kwargs) == 0 and len(args) == 0, \
             ValueError('Use of parent is exclusive')
    KExtract.__init__(self, *args, **kwargs)
    if parent is not None: self.__dict__.update(parent.__dict__)
  
  @property
  @FileCache('LDOSCAR')
  def raw_ldos(self):
    """ Raw Local density of states for given sets of positions. """
    from ldos import ldos as outer_ldos
    return outer_ldos(self, self.positions, raw=True)

  @property
  def positions(self):
    """ Positions for which to compute LDOS. """
    from numpy import array
    from quantities import angstrom
    if getattr(self.functional, 'positions', None) is None:
      return array([a.pos * self.structure.scale * angstrom for a in self.structure.atoms])
    if not hasattr(self.functional.positions, '__call__'): return self.functional.positions
    return self.funtional.positions(self.structure)

  @property
  @make_cached
  def ldos(self):
    """ Local density of states for ``positions``. """
    from numpy.linalg import det, inv
    volume = det(inv(self.structure.cell))
    return _ldosfunc(self.copy(unreduce=False).eigenvalues.flatten(), self.raw_ldos, volume)
   
  def iterfiles(self, **kwargs):
    """ Iterates through exportable files. 

        All parameters passed on to KExtract. 
        Adds LDOSCAR to export files.
    """ 
    from os.path import exists, join
    path = join(self.directory, 'LDOSCAR')
    if exists(path): yield path
    for file in KExtract.iterfiles(self, **kwargs): yield file

  @property 
  def success(self):
    """ True if successful run. """
    from os.path import join, exists
    if not exists(join(self.directory, 'LDOSCAR')): return False
    return KExtract.success.__get__(self)

class Functional(KEscan): 
  """ Functional to compute local density of states. """
  Extract = Extract
  """ Extraction object for LDOS. """
  def __init__(self, **kwargs):
    """ Initializes an LDOS functional. 
    
        :param kwargs: Any keyword argument that works for `KEscan`.
        :kwarg positions: callable which takes a structure and returns an
          array of positions where to perform ldos. Can be None, in which case,
          it defaults to the atomic positions. Must also be pickleable.
    """
    super(Functional, self).__init__(**kwargs)
    self.positions = kwargs.pop('positions', None)
    """ Callable returning positions for local density of states.

        Should be None (in which case atomic positions are used) or a
        pickleable callable which takes a structure and returns the positions
        for which to compute LDOS.
    """

  def __call__(self, *args, **kwargs):
    """ Calls KEscan, and then calls LDOS itself. 

        All parameters are passed on to escan.
    """ 
    out = super(Functional, self).__call__(*args, **kwargs)
    result = self.Extract(parent=out)
    result.ldos # computes and saves ldos.
    return result



      
