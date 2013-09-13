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

__docformat__ = "restructuredtext en"
__all__ = ['relax', 'iter_relax', 'Relax', 'RelaxExtract']
from ..tools.makeclass import makeclass, makefunc
from .functional import Functional
from .extract import Extract, MassExtract

class RelaxExtract(Extract):
  """ Extractor class for vasp relaxations. """
  class IntermediateMassExtract(MassExtract):
    """ Focuses on intermediate steps. """
    def __init__(self, *args, **kwargs):
      """ Makes sure we are using ordered dict. """
      from ..jobfolder.ordered_dict import OrderedDict
      super(RelaxExtract.IntermediateMassExtract, self).__init__(*args, **kwargs)
      self.dicttype = OrderedDict
      """ Type of dictionary to use. 
      
          Always ordered dictionary for intermediate mass extraction.
          Makes it easier to explore ongoing calculations.
      """

    def __getitem__(self, value):
      """ Gets intermediate step. 

          Adds ability to get runs as numbers: ``self.details[0]`` 
          is equivalent to ``self.details['relax/0']``
      """
      from ..error import IndexError
      if isinstance(value, int): 
        N = len(self)
        if value < 0: value += N
        if value < 0 or value >= N:
          raise IndexError( 'details: Index out of range -- no run {0}.'       \
                            .format(value) )
        value = '/relax/{0}'.format(value)
      return super(RelaxExtract.IntermediateMassExtract, self)                 \
                .__getitem__(value)

    def __iter_alljobs__(self):
      """ Goes through all directories with an OUTVAR. """
      from re import match
      from glob import iglob
      from os.path import relpath, join, exists

      results = []
      for dir in iglob(join(join(self.rootpath, 'relax'), '*/')):
        if not match('\d+', dir.split('/')[-2]): continue 
        if not exists(join(self.rootpath, join(dir, 'crystal.out'))): continue
        try: result = Extract(dir[:-1])
        except: continue
        results.append((join('/', relpath(dir[:-1], self.rootpath)), result))
      results = sorted(results, key=lambda x: int(x[0].split('/')[-1]))
      return results

    def _complete_output(self, structure):
      """ Completes output after stopped jobs. """
      result = False
      dummy = structure
      for extractor in self.itervalues():
        if hasattr(extractor, '_complete_output'):
          if extractor._complete_output(dummy):
            extractor.uncache()
            result = True
        dummy = extractor.input_crystal
      return result
        
  @property
  def details(self):
    """ Intermediate steps. """
    if '_details' not in self.__dict__:
      self.__dict__['_details'] = self.IntermediateMassExtract(self.directory)
    return self._details

  @property
  def last_step(self):
    """ Extraction object for current step. """
    if self.success: return self
    def cmp(a): return int(a[0].split('/')[-1])
    return max(self.details.items(), key=cmp)[1]
  
  def iterfiles(self, **kwargs):
    """ Iterates over input/output files. """
    from itertools import chain
    for file in chain( super(RelaxExtract, self).iterfiles(**kwargs),
                       self.details.iterfiles(**kwargs) ): yield file

  def _complete_output(self, structure):
    """ Completes output after stopped jobs. """
    result = self.details._complete_output(structure)
    try: structure = self.details[-1].input_crystal
    except: pass
    else: 
      if super(RelaxExtract, self)._complete_output(structure):
        return True
    return result

  @property
  def is_running(self):
    """ True if program is running on this functional. 
         
        A file '.pylada_is_running' is created in the output folder when it is
        set-up to run CRYSTAL_. The same file is removed when CRYSTAL_ returns
        (more specifically, when the :py:class:`pylada.process.ProgramProcess` is
        polled). Hence, this file serves as a marker of those jobs which are
        currently running.
    """
    from os.path import join, exists
    if exists(join(self.directory, '.pylada_is_running')): return True
    for value in self.details.itervalues():
      if value.is_running: return True
    return False

def iter_relax(self, structure=None, outdir=None, maxiter=30, **kwargs):
  """ Performs relaxations until convergence is reached.
  
      A final static calculation is performed at the end of the run.
  """
  from os.path import join
  from os import getcwd
  from ..error import ExternalRunFailed

  self = self.copy()
  if outdir is None: outdir = getcwd()
  self.optgeom.enabled = True
  if self.optgeom.maxcycle is None: 
    self.optgeom.maxcycle = 10
  elif self.optgeom.maxcycle < 0: 
    self.optgeom.maxcycle = 10
   
  
  iteration, restart = 0, kwargs.pop('restart', None)

  while maxiter > 0 and iteration < maxiter:
    # performs calculation
    outpath = join(join(outdir, 'relax'), str(iteration))
    for extract in self.iter( structure, outdir=outpath, restart=restart,
                              **kwargs ):
      yield extract

    # check for success
    if not extract.success: 
      raise ExternalRunFailed( 'CRYSTAL run did not succeed in {0}.'           \
                               .format(outpath) )

    # otherwise, go to next iteration
    structure = extract.crystal
    restart   = extract
    iteration += 1

    # check for convergence
    # last condition makes it easier to increase tolerance after first few run
    # at lower toldee.
    if extract.optgeom_convergence is True                                     \
       and extract.optgeom_iterations == 1                                     \
       and extract.params.toldee >= self.toldee: break

  if extract.optgeom_convergence is True                                       \
     and extract.optgeom_iterations == 1                                       \
     and extract.params.toldee >= self.toldee:
    # perform static calculation
    self.optgeom.enabled = False
    self.guessp = False
    for extract in self.iter( structure, outdir=outdir, restart=restart, 
                              **kwargs ):
      yield extract

  yield iter_relax.Extract(outdir)     

iter_relax.Extract = RelaxExtract
""" Extraction object for relaxation meta-functional. """

Relax = makeclass( 'Relax', Functional, iter_relax, None,
                   module='pylada.dftcrystal.relax',
	            	   doc='Functional form of the '                               \
                       ':py:class:`pylada.dftcrystal.relax.iter_relax` method.' )
relax = makefunc('relax', iter_relax, module='pylada.dftcrystal.relax')
