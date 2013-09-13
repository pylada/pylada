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

""" Interface module for ESCAN. """
__docformat__ = "restructuredtext en"
__all__ = [ "Extract", 'MassExtract', "bandgap", "extract_bg", 'ldos',
            "Escan", "folded_spectrum", "all_electron", "soH", 'BandGap',
            "nonlocalH", "localH", "AtomicPotential", "extract_bg", 'KEscan', 'KPoints', 
            'KGrid', 'ReducedKGrid', 'ReducedKDensity', 'soH', 'nonlocalH', 'localH', 
            'folded_spectrum', 'all_electron', 'read_input', 'exec_input', 'KExtract',  
            'majority_representation', 'BPoints', 'ReducedBPoints', 'plot_bands', 'plot_alloybands',
            'fftmesh', 'EMassFunctional', 'EMassExtract', 'InnerBPoints', 'ReducedInnerBPoints' ]

from _bandstructure import plot_bands, BPoints, ReducedBPoints,\
                           plot_alloybands, InnerBPoints, ReducedInnerBPoints
from _bandgap import bandgap, extract as extract_bg, Functional as BandGap
from emass import Functional as EMassFunctional, Extract as EMassExtract
from _extract import Extract
from _massextract import MassExtract
from _potential import soH, nonlocalH, localH, AtomicPotential
from _methods import majority_representation
from functional import Escan, folded_spectrum, all_electron
from kescan import KEscan, Extract as KExtract
from kpoints import KGrid, ReducedKGrid, KPoints, ReducedKDensity
import ldos
import fftmesh


def exec_input(script, namespace = None):
  """ Executes an input script including namespace for escan/vff. """ 
  from ..opt import exec_input as opt_exec_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace is not None: dictionary.update(namespace)
  return opt_exec_input(script, dictionary)

def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..opt import read_input as opt_read_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace is not None: dictionary.update(namespace)
  return opt_read_input(filepath, dictionary)

def effective_mass_tensor(escan, structure, outdir=None, comm=None, attenuate=None, degeneracy=None, **kwargs):
  """ Computes effective mass tensor using dipoles and f-sum rule. 
  
      :Parameters:
        escan : escan.Escan
          Straight off escan functional. Calculations must be full diagonalization.
        structure : pylada.crystal.Structure
          Crystal structure for which to perform calculations.
        outdir : str or None
          Output directory. Defaults to current directory.
        comm : pylada.mpi.communicator.
          MPI communicator with which to perform calculation. 

      :return: a numpy array with units of an electronic mass, with the 
               dimensions of the input directions, except for the last axis,
               which is of the same length as the number of computed
               eigenvalues.
  """
  if kwargs.get('eref', escan.eref) != None: 
    raise ValueError("f-sum rule computation of effective mass only done with direct diagonalization.")

  # compute escan wavefunctions.
  result = escan(structure, outdir, comm, **kwargs)

  # compute effective mass tensor.
  ekwargs = {}
  if attenuate != None: ekwargs['attenuate'] = attenuate
  if degeneracy != None: ekwargs['degeneracy'] = degeneracy
  return result.effective_mass_tensor(**ekwargs)

effective_mass_tensor.Extract = Extract
""" Extraction object for f-sum rule effective masses. """
