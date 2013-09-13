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
from ..tools import stateless, assign_attributes
from ..tools.input import AttrBlock
from .extract import Extract as ExtractBase

class Functional(AttrBlock):
  """ Wrapper around the GULP program. """
  Extract = ExtractBase
  """ Extraction class. """
  __ui_name__ = 'functional'
  """ Name used in user-friendly representation """

  def __init__(self, program=None, copy=None, **kwargs):
    """ Creates the GULP wrapper. """
    from ..tools.input.keywords import BoolKeyword
    from .keywords import TwoBody, OptiKeyword, Optimize, Conv, Conp, Species, \
                          Springs
    super(Functional, self).__init__()
    
    self.program = program
    """ Path to crystal program.

        If this attribute is None, then :py:data:`~pylada.gulp_program` is used.
    """ 
    self.opti = Optimize()
    """ If True, then performs structural optimization. """
    self.cellonly = OptiKeyword()
    """ If True, optimization will run on cell degrees of freedom only.
    
        This keywords prints to the GULP input only if :py:attr:`opti` is on.
    """
    self.isotropic = OptiKeyword()
    """ If True, optimization will of cell is isotropic (no shear?).
    
        This keywords prints to the GULP input only if :py:attr:`opti` is on.
    """
    self.conv = Conv()
    """ If True, then performs constant volume optimization.
    
        This keywords prints to the GULP input only if :py:attr:`opti` is on.
    """
    self.conp = Conp()
    """ If True, then performs constant pressure optimization.
    
        This keywords prints to the GULP input only if :py:attr:`opti` is on.
    """
    self.shell = OptiKeyword()
    """ If True, optimizes only shells internal degrees of freedom. 
    
        This keywords prints to the GULP input only if :py:attr:`opti` is on.
    """
    self.qeq = BoolKeyword()
    """ Charge equalization scheme or Rappe and Goddard III. """
    self.morse = TwoBody()
    """ Holds parameters of the two-body Morse interaction. """
    self.buckingham = TwoBody()
    """ Holds parameters of the two-body Buckingham interaction. """
    self.species = Species()
    """ Atomic parameters. """
    self.springs = Springs()
    """ Springs between cores and shells. """

  def input_string(self, structure=None, **kwargs):
    """ Returns string string with input. """
    from ..error import ValueError
    from ..crystal.write import gulp
    # first gets structure.
    if structure is not None: 
      result = gulp(structure).rstrip()
      if len(result): result += '\n\n'
    else: result = ""
    # then computes map.
    map = self.output_map(structure=structure, gulp=self, **kwargs)
    if map is None: raise ValueError('No interactions enabled.')

    header = ""
    for key, value in map:
      if value is None or value is True                                        \
         or (isinstance(value, str) and len(value) == 0)                       \
         or (isinstance(value, str) and value == 'True'):
         header += str(key) + " "
      elif value is not False and value != 'False':
        result += '\n{0}\n{1}'.format(key, value)
    return header + '\n' + result

  def bringup(self, structure, outdir, workdir=None, **kwargs):
    """ Prepares files for gulp. """
    from os.path import join
    from ..tools import create_directory, prep_symlink, add_pyladarunning_marker,\
                        add_section_to_file
  
    create_directory(outdir)
    create_directory(workdir)

    input = self.input_string(structure, **kwargs) 
    add_section_to_file(outdir, 'gulp.out', 'input file', input, False)
    add_section_to_file(outdir, 'gulp.out', 'functional', repr(self))
    with open(join(outdir, 'gulp.in'), 'w') as file: file.write(input)
    with open(join(workdir, 'gulp.in'), 'w') as file: file.write(input)

    prep_symlink(outdir, workdir, 'gulp.out')
    prep_symlink(outdir, workdir, 'gulp.err')
    prep_symlink(outdir, workdir)

    add_pyladarunning_marker(outdir)

  def bringdown(self, structure, workdir, outdir):
    """ Cleansup after functional. """
    from os.path import join, exists
    from ..tools import remove_pyladarunning_marker, add_section_to_file,        \
                        remove_workdir_link

    remove_pyladarunning_marker(outdir)

    if exists(join(outdir, 'gulp.err')): 
      with open(join(outdir, 'gulp.err'), 'r') as file: string = file.read()
      add_section_to_file(outdir, 'gulp.out', 'standard error', string)

    if ExtractBase(outdir).success: remove_workdir_link(outdir)

  @stateless
  @assign_attributes(ignore=['overwrite', 'comm', 'workdir'])
  def iter( self, structure, outdir=None, workdir=None, comm=None,
            overwrite=False, **kwargs ):
    """ Yields processes with which to call GULP. """
    from os import getcwd
    from ..process.program import ProgramProcess
    from ..misc import Changedir, RelativePath
    from ..tools import OnFinish
    from .. import gulp_program
    # check for pre-existing and successfull run.
    if not overwrite:
      extract = ExtractBase(outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return

    if outdir == None: outdir = getcwd()
    if workdir == None: workdir = self.guess_workdir(outdir)

    outdir = RelativePath(outdir).path
    workdir = RelativePath(workdir).path
    with Changedir(workdir) as tmpdir: 

      # writes/copies files before launching.
      self.bringup(structure, outdir, workdir)
      dompi = comm is not None and comm['n'] > 1

      # figure out the program to launch.
      program = self.program if self.program is not None else gulp_program
      if hasattr(program, '__call__'):
        program = program(self, structure, comm=comm)

      # now creates the process, with a callback when finished.
      onfinish = OnFinish(self, structure, workdir, outdir)
      yield ProgramProcess( program, outdir=workdir, onfinish=onfinish,
                            stdout=('gulp.out', 'a'), stderr='gulp.err',
                            stdin='gulp.in', dompi=dompi )
    # yields final extraction object.
    yield ExtractBase(outdir)

  def __call__( self, structure, outdir=None, workdir=None, comm=None,         \
                overwrite=False, **kwargs):
    from ..error import ExternalRunFailed
    for program in self.iter( structure, outdir=outdir, workdir=workdir,
                              comm=comm, overwrite=overwrite, **kwargs ):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False): continue
      # Or may fail return a failed run.
      if not hasattr(program, 'start'): return program
      # otherwise, it should yield a Program tuple to execute.
      program.start(comm)
      program.wait()
    # Last yield should be an extraction object.
    if not program.success:
      raise ExternalRunFailed("GULP failed to execute correctly.")
    return program
  __call__.__doc__ = iter.__doc__

  def __repr__(self, defaults=True, name=None):
    """ Representation of this instance. """
    from ..tools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def guess_workdir(self, outdir):
    """ Tries and guess working directory. """
    from ..misc import mkdtemp
    from .. import crystal_inplace
    return outdir if crystal_inplace else mkdtemp(prefix='gulp') 
