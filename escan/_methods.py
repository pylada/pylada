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

""" Subpackage containing dipole matrix-element implementations. """
__docformat__ = "restructuredtext en"


def majority_representation(extractor, multicell, tolerance=1e-12): 
  """ Majority Represention for all unmapped kpoint. 

      :Parameters:
        extractor : `KExtract`
           Result of a `KEscan` calculation. The ``kpoints`` attribute from
           ``extractor.functional`` should be a class of *reduced* kpoints, such
           as `ReducedKGrid`, `ReducedKDensity`, or `ReducedBPoints`.
        multicell : numpy 3x3
           A matrix (M) linking the unit-cell (U) to the super cell (S), ``S = U
           * M``.  From this, the unit-cell is determined, whatever deformation
           may have happened, as engendered by VFF.
       
      :return: An  list of lists of 2-tuples. The first dimension corresponds
        to unreduced kpoints.  The 2-tuples consiste  of the eigenvalue for that
        kpoint and band, and of its majority representation. The inner list
        duplicate entries with respect to eigenvalues. It is meant to be used
        with matplotlib's hexbin.

      The majority representation is a way to get to the band-structure of an
      alloy despite the lack of translational symmetries. It was introduced by
      L.L.Wang *et al.* in  [LLW]_, and then further refined in [VP]_.

      .. [LLW] Lin Wang Wang, Laurent Bellaiche, Suhai Wei, and Alex Zunger,
         Phys. Rev. Lett. **80**, 4725 (1998). 

      .. [VP] Voicu Popescu and Alex Zunger, Phys. Rev. Lett. **104**, 236403
         (2010).
  """
  from operator import itemgetter
  from numpy import dot, pi, array
  from numpy.linalg import inv
  from quantities import angstrom
  from ..crystal import is_on_lattice, to_voronoi
  from ..physics import a0

  assert extractor.success,\
         ValueError("extractor is not the return of a *successfull* KEscan calculation.")
  assert len(extractor) > 0, ValueError("extractor does not contain escan runs.")
  istr = extractor.input_structure
  ostr = extractor.structure
  unitcell = dot(ostr.cell, inv(multicell)) * ostr.scale * angstrom
  unitkcell = inv(unitcell.T) * 2e0 * pi / angstrom
  invcell = inv(ostr.cell.T)
  
  # use the functional's kpoint object so we have access to all its goodies.
  kpoints = extractor.functional.kpoints

  # list of results.
  results = []

  unreduced = array([k[1] for k in kpoints.unreduced(istr, ostr)])
  mapping = [i for i in kpoints.mapping(istr, ostr)]
  scale = float(1e0/(ostr.scale  * angstrom).rescale(a0))
  N = float(len(ostr.atoms))
  for k0, extract in zip(unreduced, extractor):
    # actual calculation for that kpoint.
    # reciprocal vectors of supercell.  
    Gs = extract.gvectors        
    # Compute kpoint folded into voronoi cell.
    K = to_voronoi(k0, invcell)
    # links reduced to unreduced: folding vector.
    Gfold = (K - k0) * scale
    # mapping for reciprocal vector of supercell which are reciprocal vectors
    # of unitcell *with* folding.
    onlatop = is_on_lattice(Gs, unitkcell.rescale(Gs.units), Gfold)
    # loop over bands and sum 
    current = [[w.eigenvalue, w.expectation_value(onlatop).real] for w in extract.gwfns]
    current = sorted(current, key=itemgetter(0))
    intermediate = [current[0]]
    for eig, val in current[1:]: 
      if abs(eig-intermediate[-1][0]) < tolerance: intermediate[-1][1] += val
      else: intermediate.append([eig, val])

    results.append([ (e, v) for e, v in intermediate if abs(v) * N > tolerance])

  return results
