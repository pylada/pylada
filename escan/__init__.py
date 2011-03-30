""" Interface module for ESCAN. """
__docformat__ = "restructuredtext en"
__all__ = [ "Extract", 'MassExtract', "bandgap", "extract_bg", 'ldos',
            "Escan", "folded_spectrum", "all_electron", "soH", 
            "nonlocalH", "localH", "AtomicPotential", "extract_bg", 'KEscan', 'KPoints', 
            'KGrid', 'ReducedKGrid', 'ReducedKDensity', 'soH', 'nonlocalH', 'localH', 
            'folded_spectrum', 'all_electron', 'read_input', 'exec_input', 'KExtract',  
            'majority_representation', 'BPoints', 'ReducedBPoints', 'plot_bands', 'plot_alloybands']

from ..opt import __load_escan_in_global_namespace__
from _bandstructure import plot_bands, BPoints, ReducedBPoints, plot_alloybands
from _bandgap import bandgap, extract as extract_bg
from _extract import Extract
from _massextract import MassExtract
from _potential import soH, nonlocalH, localH, AtomicPotential
from _methods import majority_representation
from functional import Escan, folded_spectrum, all_electron
from kescan import KEscan, Extract as KExtract
from kpoints import KGrid, ReducedKGrid, KPoints, ReducedKDensity
import ldos


def exec_input(script, namespace = None):
  """ Executes an input script including namespace for escan/vff. """ 
  from ..opt import exec_input as opt_exec_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace != None: dictionary.update(namespace)
  return opt_exec_input(script, dictionary)

def read_input(filepath = "input.py", namespace = None):
  """ Reads an input file including namespace for escan/vff. """ 
  from ..opt import read_input as opt_read_input
  from .. import vff

  dictionary = {}
  for key in vff.__all__: dictionary[key] = getattr(vff, key)
  for key in __all__: dictionary[key] = globals()[key]
  if namespace != None: dictionary.update(namespace)
  return opt_read_input(filepath, dictionary)
