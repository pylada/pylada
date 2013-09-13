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

""" Subpackage grouping single-electron properties. """
__docformat__ = "restructuredtext en"
__all__ = ['Properties']
from ..input import AttrBlock
from .extract import Extract

class Properties(AttrBlock):
  """ Wrap single-electron property calculations. """
  Extract = Extract
  """ Extraction class. """
  __ui_name__ = 'properties'
  """ Name used in user-friendly representation """
  def __init__(self, input=None, program=None):
    from .keywords import NewK, Band, Rdfmwf, Fmwf
    from ..input import BoolKeyword
    super(Properties, self).__init__()
    self.rdfmwf = Rdfmwf()
    """ Whether to create binary wavefunction file. 
    
        - If None, Pylada figures out which of  'crystal.f98' and 'crystal.f9'
          file is the latest in input and output directories. Then figures out
          whether to include an RDFMWF tag or not.
        - If True, Pylada looks only for a 'crystal.f98' in the input and/or
          output directories, using the latest. If no file exists, then it is
          an error.
        - If False, Pylada checks that 'crystal.f9' exists in the input and/or
          output directories, using the latest. If no file exists, then it is
          an error.

        If no error occurs, the relevant files are copied to the working
        directory.

        .. note:: 

           Only one of :py:attr:`rdfmwf` and :py:attr:`fmwf` can be True at any
           one time.
    """
    self.newk = NewK()
    """ Performs diagonalization on new k-point mesh. 

        Starting from the hamiltonian defined in input, performs a
        diagonalization on a new k-point mesh.
        The input is similar to :py:attr:`~pylada.functional.shrink`.
        Additionaly, it is possible to specify printing options:

        >>> properties.newk.printing[66] = -5

        Following the same input scheme as
        :py:attr:`~pylada.functional.setprint`.

        It is possible to set the k-point mesh directly, much as for
        :py:attr:`~pylada.functional.shrink`:

        >>> properties.newk = 5, None

        However, both the ``recompute_fermi`` and the printing options
        ``printing`` must be set by direct access:

        >>> properties.newk.printing[5] = 0
        >>> properties.recompute_fermi = True
    """
    self.nosymada = BoolKeyword()
    """ No symmetry adapted bloch function. """
    self.band = Band()
    """ Performs calculation on a path of k-points.
    
        This attribute can be parameterized in a variety of way.

        - ``properties.band.title``: accepts a string which will be the title
          of the band-structure. This is mostly to reproduce original CRYSTAL_
          input.
        - ``properties.band.nsub``: accepts an integer which defines the total
          number of points along the path.
        - ``properties.band.iss``: see CRYSTAL_'s user guide. Not really
          relevant in python.
        - ``properties.band.minband``: Index of the minimum band for which to
          output results. Best leave this to zero since bands starting from
          zero are computed anyways. Hold-out from CRYSTAL_.
        - ``properties.band.maxband``: Index of the maximum band for which to
          ouput results.

        The path is held the ``lines`` attributes. It can be given as:

        >>> properties.band = [startA, endA], [startB, endB], ...
    """
    self.fmwf = Fmwf()
    """ Whether to output formatted wavefunctions. 

        .. note:: 

           Only one of :py:attr:`rdfmwf` and :py:attr:`fmwf` can be True at any
           one time.
    """
    self.program = program
    """ CRYSTAL_'s properties program. 

        If None, defaults to :py:data:`pylada.properties_program`.
    """ 
    self.input = input

  def output_map(self, **kwargs):
    """ Returns map of crystal input. """
    from ...tools.input import Tree
    from ...misc import RelativePath
    if 'properties' not in kwargs: kwargs['properties'] = self
    if 'input' not in kwargs: kwargs['input'] = self.input
    if 'outdir' not in kwargs: kwargs['outdir'] = None
    kwargs['outdir'] = RelativePath(kwargs['outdir']).path 
    if 'workdir' not in kwargs: kwargs['workdir'] = None
    kwargs['workdir'] = RelativePath(kwargs['workdir']).path 
    root = Tree()
    # First add RDFWF
    AttrBlock._output_map(root, 'rdfmwf', self._input['rdfmwf'], **kwargs)
    # Then go on to symmetry adapted thingie
    AttrBlock._output_map(root, 'nosymada', self._input['nosymada'], **kwargs)
    # Move on to newk
    AttrBlock._output_map(root, 'newk', self._input['newk'], **kwargs)
    # Do other preliminaries.
    for key in ['pato', 'pban', 'pgeomw', 'pdide', 'pscf']:
      if key in self._input: 
        AttrBlock._output_map(root, key, self._input[key], **kwargs)
      elif key.upper() in self._input:
        AttrBlock._output_map(root, key.upper(), self._input[key], **kwargs)
    # Now do all others.
    prelims = set([ 'rdfmwf', 'nosymada', 'newk', 'pato', 
                    'pban', 'pgeomw', 'pdide', 'pscf' ]) 
    for key, value in self._input.iteritems():
      if key not in prelims and key.upper() not in prelims: 
        AttrBlock._output_map(root, key, value, **kwargs)
    return root
  
  def print_input(self, **kwargs):
    """ Prints input to string. """
    from ..input import print_input
    map = self.output_map(**kwargs)
    if len(map) == 0: return ""
    # Otherwise, everything is standard.
    return print_input(map).rstrip() + '\nEND\n'

  @property
  def input(self):
    """ Input calculation from which to start. 

        - If None, then an extraction object is created in the current
          directory.
        - If a string, then it is assumed to the path to the directory
          containing the input calculations.
        - If an extraction object, then that extraction object should point to
          the self-consistent calculations.
    """
    return self._input_calc
  @input.setter
  def input(self, value):
    from ..extract import Extract as CrystalExtract
    if value is None: self._input_calc = CrystalExtract()
    elif isinstance(value, str): self._input_calc = CrystalExtract(value)
    else: self._input_calc = value

  def nbelectrons(self, input=None):
    if input is None: input = self.input
    species = [u.type for u in input.structure]
    result = 0
    for specie in set(species):
      result += sum(u.charge for u in input.functional.basis[specie])          \
                * species.count(specie)
    return result

  def guess_workdir(self, outdir):
    """ Tries and guess working directory. """
    from ...misc import mkdtemp
    from ... import crystal_inplace
    return outdir if crystal_inplace else mkdtemp(prefix='crystalprop') 

  def iter( self, input=None, outdir=None, workdir=None, overwrite=False,
            **kwargs ):
    """ Performs a single-electron property calculation. """
    from os import getcwd
    from os.path import abspath
    from copy import deepcopy
    from ...process.program import ProgramProcess
    from ...misc import Changedir
    from ... import properties_program
    from ..functional import Functional
    from ...error import AttributeError, input as InputError

    # Performs deepcopy of self.
    this = deepcopy(self)
    # copy attributes.
    if input is not None: this.input = input
    for key, value in kwargs.iteritems():
      if key in ['comm', 'structure']: continue
      if hasattr(self, key): setattr(self, key, value)
      else: raise AttributeError('Unknown attribute {0}'.format(key))
      
    # check for pre-existing and successfull run.
    if not overwrite:
      extract = this.Extract(outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return
    
    if outdir == None: outdir = getcwd()
    # Makes sure that we have actual work to do.
    # we check prior to checking the workdir since that might create it.
    map = self.print_input(workdir=outdir, outdir=outdir, filework=False)
    if len(map) == 0:
      raise InputError( 'Nothing to do. Properties not requested '           \
                        'to compute anything.' )
    # now make sure guessdir is well defined.
    if workdir is None: workdir = this.guess_workdir(outdir)

    outdir = abspath(outdir)
    workdir = abspath(workdir)
    with Changedir(workdir) as tmpdir: 
      # writes/copies files before launching.
      this.bringup(outdir, workdir)

      # figure out the program to launch.
      program = this.program if this.program is not None else properties_program
      if hasattr(program, '__call__'): program = program(this)

      # now creates the process, with a callback when finished.
      onfinish = Functional.OnFinish(this, workdir, outdir)
      yield ProgramProcess( program, outdir=workdir, onfinish=onfinish,
                            stdout=('prop.out', 'a'), stderr='prop.out', 
                            stdin='prop.d12', dompi=False )
    # yields final extraction object.
    yield Extract(outdir)

  def bringup(self, outdir, workdir):
    """ Sets up call to program. """
    from os.path import join, abspath, samefile, lexists
    from os import symlink, remove
    from ...misc import copyfile, Changedir
    from ... import CRYSTAL_propnames as filenames
    with Changedir(workdir) as cwd:
      # first copies file from current working directory
      for key, value in filenames.iteritems():
        copyfile( join(workdir, value.format('prop')), key, nocopyempty=True,
                  symlink=False, nothrow="never" )
      for key, value in filenames.iteritems():
        copyfile( join(self.input.directory, value.format('prop')), 
                  key, nocopyempty=True, symlink=False, 
                  nothrow="never" )

      # then creates input file.
      string = self.print_input(workdir=workdir, outdir=outdir, filework=True)
      string = string.rstrip() + '\n'
      with open('prop.d12', 'w') as file: file.write(string)
      header = ''.join(['#']*20)
      with open('prop.out', 'w') as file:
        file.write('{0} {1} {0}\n'.format(header, 'INPUT FILE'))
        file.write(string)
        file.write('{0} END {1} {0}\n'.format(header, 'INPUT FILE'))
        file.write('\n{0} {1} {0}\n'.format(header, 'FUNCTIONAL'))
        file.write(self.__repr__(defaults=False))
        file.write('\n{0} END {1} {0}\n'.format(header, 'FUNCTIONAL'))

    with Changedir(outdir) as cwd: pass
    if not samefile(outdir, workdir):
      # Creates symlink to make sure we keep working directory.
      with Changedir(outdir) as cwd:
        with open('prop.d12', 'w') as file: file.write(string)
        with open('prop.out', 'w') as file: pass
        with open('prop.err', 'w') as file: pass
        # creates symlink files.
        for filename in ['prop.err', 'prop.out']:
          if lexists(join(workdir, filename)):
            try: remove( join(workdir, filename) )
            except: pass
          symlink(abspath(filename), abspath(join(workdir, filename)))
            
        if lexists('workdir'): 
          try: remove('workdir')
          except: pass
        try: symlink(workdir, 'workdir')
        except: pass

    # creates a file in the directory, to say we are going to work here
    with open(join(outdir, '.pylada_is_running'), 'w') as file: pass

  def bringdown(self, workdir, outdir):
    """ Copies files back to output directory. 
    
        Cats input intO output. Removes workdir if different from outdir
        **and** run was successfull.
    """
    from os import remove
    from os.path import join, exists, samefile
    from shutil import rmtree
    from ...misc import copyfile, Changedir
    from ... import CRYSTAL_propnames as propnames,                            \
                    CRYSTAL_filenames as crysnames

    with Changedir(outdir) as cwd:
      for key, value in propnames.iteritems():
        copyfile( join(workdir, key), value.format('prop'),
                  nocopyempty=True, symlink=False, nothrow="never" )
      if self.fmwf is True:
        copyfile( join(workdir, 'fort.98'),
                  crysnames['fort.98'].format('crystal'),
                  nocopyempty=True, symlink=False, nothrow="never" )
      if self.rdfmwf is True:
        copyfile( join(workdir, 'fort.9'),
                  crysnames['fort.9'].format('crystal'),
                  nocopyempty=True, symlink=False, nothrow="never" )

      # remove 'is running' file marker.
      if exists('.pylada_is_running'):
        try: remove('.pylada_is_running')
        except: pass
    
    if Extract(outdir).success:
      if exists(workdir) and not samefile(workdir, outdir):
        try: rmtree(workdir)
        except: pass
      try: remove(join(outdir, 'workdir'))
      except: pass

  def __call__( self, input=None, outdir=None, workdir=None, overwrite=False,
                **kwargs):
    for program in self.iter( input=input, outdir=outdir, workdir=workdir,
                              overwrite=overwrite, **kwargs ):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False): continue
      # Or may fail return a failed run.
      if not hasattr(program, 'start'): return program
      # otherwise, it should yield a Program tuple to execute.
      program.start()
      program.wait()
    # Last yield should be an extraction object.
    if not program.success:
      raise RuntimeError("CRYSTAL failed to execute correctly.")
    return program
  __call__.__doc__ = iter.__doc__

  def __repr__(self, defaults=True, name=None): 
    """ Returns representation of this instance """
    from ...tools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)
