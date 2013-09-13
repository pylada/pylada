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

""" Contains evaluators for ESCAN properties """
__docformat__ = "restructuredtext en"
from ...opt.decorators import count_calls

__all__ = ['Bandgap', 'Dipole', 'EffectiveMass']

  
class Bandgap(object):
  """ An evaluator function for direct bandgaps. """
  def __init__(self, converter, functional, outdir = None, references = None, \
               keep_only_last = True, sizedependent = -1, **kwargs):
    """ Initializes the bandgap object. 

        :Parameters: 
          converter : `pylada.ga.escan.elemental.Converter`
            is a functor which converts between a bitstring and a
            `pylada.crystal.Structure`, and vice versa.
          functional 
            some escan functional or meta-functional.
          outdir : str or None
            Directory where to perform calculations. Defaults to None. 
          references
	          Reference energies for band-gap calculations. See
            `pylada.escan.bandgap`. Defaults to None.
          keep_only_last : boolean
	          If true, only the last calculation is kept. Limits the space
            required for ESCAN optimizations.
          sizedependent : float
            If strictly postive, then number of states is the product of this
            number and the number of bits in the bitstring.
          kwargs
            Other keyword arguments to be passed to the functional.
    """
    from copy import deepcopy
    from ...opt import RelativeDirectory

    super(Bandgap, self).__init__()

    self.converter = converter 
    """ Conversion functor between structures and bitstrings. """
    self.escan = functional
    """ Escan functional or meta-functional. """
    self.keep_only_last = keep_only_last
    """ Whethere to keep only last calculations or all. """
    self.kwargs = deepcopy(kwargs)
    """ Additional arguments to be passed to the bandgap functional. """
    self.sizedependent = sizedependent
    """ Floating points to compute a size-dependent number of states. 
 
        If strictly postive, then number of states is the sum of escan.nbstates
        with  the product between this number and the size of the bitstring.
    """


    self._outdir = RelativeDirectory(path=outdir)
    """ Location of output directory. """
    self._lastcalcdir = None
    """ Location of last calculation. """

  @property 
  def outdir(self):
    """ Current workding directory. """
    return self._outdir.path
  @outdir.setter 
  def outdir(self, value): self._outdir.path = value

  def __len__(self):
    """ Returns length of bitstring. """
    if hasattr(self.converter, '__len__'): return len(self.converter)
    raise AttributeError("{0} has not len attribute.".format(self.__class__.__name__))

  @count_calls('nbcalc', 0)
  def run(self, indiv, outdir = None, comm = None, **kwargs):
    """ Computes bandgap of an individual. 
    
        :kwarg indiv: Individual to compute. Will be converted to a
          `pylada.crystal.Structure` using `converter`.
        :kwarg outdir: Output directory. 
        :kwarg comm: `mpi.Communicator`.
        :kwarg kwarg: converter, escan, references, outdir  can be overwridden
            on call. This will not affect calls further down the line. Other
            arguments are passed on to the pylada.escan.bandgap` function on call.

        :return: an extractor to the bandgap. 

        The epitaxial energy is stored in indiv.epi_energy
        The eigenvalues are stored in indiv.eigenvalues
        The VBM and CBM are stored in indiv.bands
    """
    from shutil import rmtree
    from os.path import join, exists
    from copy import deepcopy
    from ...mpi import Communicator

    # takes into account input arguments.
    converter  = kwargs.pop( "converter",  self.converter)
    escan      = kwargs.pop("functional", self.escan)
    comm       = Communicator(comm)
    if outdir is None: outdir = self.outdir
    elif outdir[0] != '/': outdir = join(self.outdir, outdir)
    outdir = join(outdir, str(self.nbcalc))
 
    # creates a crystal structure (phenotype) from the genotype.
    structure = converter(indiv.genes)
    # creates argument dictonary to pass on to calculations.
    dictionary = deepcopy(self.kwargs)
    dictionary.update(kwargs) 
    dictionary["comm"]       = comm 
    dictionary["outdir"]     = outdir
    if "overwrite" not in dictionary: dictionary["overwrite"]  = True
    if getattr(self, "sizedependent", -1) > 0.0:
      dictionary["nbstates"] = float(len(indiv.genes)) * self.sizedependent + self.escan.nbstates
      dictionary["nbstates"] = 2*(int(dictionary["nbstates"]+1.5)//2)
    # performs calculation.
    out = escan(structure, **dictionary)

    # saves some stuff
    indiv.epi_energy = out.energy
    indiv.stress = out.stress.copy()
    indiv.bandgap = out.bandgap
    indiv.vbm = out.vbm
    indiv.cbm = out.cbm

    comm.barrier()
    if self.keep_only_last and comm.is_root and self._lastcalcdir is not None:
      if exists(self._lastcalcdir): rmtree(self._lastcalcdir)
    self._lastcalcdir = outdir
    comm.barrier()
    
    # returns extractor
    return out

  def __call__(self, *args, **kwargs):
    """ Computes and returns bandgap. 
    
        see self.run(...) for more details.
    """
    return self.run(*args, **kwargs).bandgap

  def __getstate__(self):
    from pickle import dumps
    d = self.__dict__.copy()
    try:  dumps(d.get("references", None))
    except TypeError: 
      raise RuntimeError("Cannot pickle references in Bandgap evaluator.")
    return d
  def __setstate__(self, arg):
    self.__dict__.update(arg)

  def __repr__(self): 
    """ Returns representation of evaluator. """
    from operator import itemgetter
    max_length, string, _string, values = 0, '', '', {}
    for key, value in sorted(self.__dict__.items(), key=itemgetter(0)):
      if key[0] == '_': continue
      if key == 'converter': continue
      if key == 'escan': continue
      if key == 'outdir': continue
      try: r = repr(value).rstrip().lstrip()
      except: continue
      else: r = r.replace('{', '{{').replace('}', '}}')
      if r[0] == '<' or r[-1] == '>': continue
      max_length = max(max_length, len('{0}'.format(key)))
      string += 'evaluator.{{{0}: <{{_mxlgth_repr_}}}} = {1}\n'.format(key, r)
      values[key] = key
    # create format string for private data members.
    for key, value in sorted(self.__dict__.items(), key=itemgetter(0)):
      if key[0] != '_': continue
      if key == '_outdir': continue
      try: r = repr(value).rstrip().lstrip()
      except: continue
      else: r = r.replace('{', '{{').replace('}', '}}')
      if r[0] == '<' or r[-1] == '>': continue
      max_length = max(max_length, len('{0}'.format(key)))
      _string += 'evaluator.{{{0}: <{{_mxlgth_repr_}}}} = {1}\n'.format(key, r)
      values[key] = key
    values['_mxlgth_repr_'] = max_length
    string += "evaluator.{{outdir: <{{_mxlgth_repr_}}}} = {1}\n"\
              .format('outdir', self._outdir.repr())
    values['outdir'] = 'outdir'

    result = "from {0.__class__.__module__} import {0.__class__.__name__}\n"\
             "{1}\n"\
             "{2}\n"\
             "evaluator = {0.__class__.__name__}(converter, escan_functional)\n"\
             .format( self, repr(self.escan).replace("functional", "escan_functional"),
                      repr(self.converter) )
    result += string.format(**values)
    result += _string.format(**values)
    return result


class Dipole(Bandgap):
  """ Evaluates the oscillator strength.

      On top of those quantities saved by base class BandgapEvaluator,
      this class stores the dipole elements in indiv.dipoles.
  """
  def __init__(self, *args, **kwargs): 
    """ Initializes the dipole element evaluator.  """
    self.degeneracy = kwargs.pop("degeneracy", 1e-3)
    """ Degeneracy parameter for oscillator strength. """
    super(Dipole, self).__init__(*args, **kwargs)

  def __call__(self, indiv, outdir=None, comm=None, **kwargs):
    """ Computes the oscillator strength. """
    out = super(Dipole, self).run(indiv, outdir, comm, **kwargs)
    degeneracy = getattr(self, 'degeneracy', 1e-3)
    indiv.dipoles = out.dipole(degeneracy=degeneracy)
    indiv.oscillator_strength, indiv.osc_nbstates\
      = out.oscillator_strength(degeneracy=degeneracy)
    comm.barrier() 
    return indiv.oscillator_strength / float(len(indiv.genes))\
           if getattr(self, "sizedependent", -1)  > 0e0 \
           else indiv.oscillator_strength


class EffectiveMass(Bandgap):
  """ Evaluates effective mass. """
  def __init__( self, n=0, **kwargs ):
    """ Computes effective mass for a given direction. 

        For other keyword arguments, see `Bandgap`.
    
        :Parameters:
          n : int or callable
            Index of the band for which to compute effective mass.
            If callable, takes the extraction object on input and should return
            a number.
    """
    self.n = n
    """ Index of the band for which to compute effective mass.

        Can also be a callable taking the extraction object on input and
        returning a number.
    """
    super(EffectiveMass, self).__init__(**kwargs)

  def __call__(self, indiv, outdir=None, comm=None, **kwargs):
    """ Computes electronic effective mass. """
    from numpy import average
    out = super(EffectiveMass, self).run(indiv, outdir, comm, **kwargs)
    indiv.eigenvalues = out.eigenvalues
    indiv.masses = out.mass
    if hasattr(self.n, "__call__"): indiv.mass = self.n(out)
    elif out.mass.ndim == 2: indiv.mass = average(out.mass[:,self.n])
    else: indiv.mass = out.mass[self.n]
    return indiv.mass * 10e0 if indiv.mass < 0e0 else indiv.mass
