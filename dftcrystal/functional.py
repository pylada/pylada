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
from .extract import Extract as ExtractBase

class Functional(object):
  """ Wrapper for the CRYSTAL program. 
  
      This object provides a pythonic interface to the CRYSTAL_ program. It is
      modelled loosely after CRYSTAL_'s input:

        - The OPTGEOM keyword of the first code block can be accessed through
          the :py:class:`optgeom <pylada.dftcrystal.optgeom.OptGeom>`
          :py:attr:`attribute <optgeom>`:
          
          .. code-block:: python

            # enable geometry optimization
            functional.optgeom.enabled = True
            # change some of the keywords
            functional.optgeom.maxcycle = True

          Geometry optimization must be explicitly enabled, as is done in the
          first line above. 
        - The interface to the second block of input (basis-functions) can be
          accessed through the :py:class:`basis
          <pylada.dftcrystal.basis.BasisSet>` :py:attr:`attribute <basis>`. It
          allows to set the basis set itself, as well as keywords specific to
          the basis set:

          .. code-block:: python
            
            functional.basis['H'] = [Shell('s', a0=(15.0, 1.0), a1=(8.0, 1.0)), ...]
            functional.basis.ghosts = [1, 4, 6], False

        - The third input block (Hamiltonian and miscellaneous) can be accessed
          directly through the functional, or, alternatively, *via* the
          :py:class:`scf <pylada.dftcrystal.electronic.Electronic>`
          :py:attr:`attribute <scf>`.
        
          .. code-block:: python
        
            functional.scf.dft.b3lyp = True
            functional.dft.b3lyp = True
        
            functional.tolinteg = [8] * 4 + [14]
        
          The first two lines are exactly equivalent. The third line could also
          be written as ``functional.scf.tolinteg = ...``.

  """
  Extract = ExtractBase
  """ Extraction class. """
  __ui_name__ = 'functional'
  """ Name used in user-friendly representation """
  def __init__(self, program=None, copy=None, **kwargs):
    """ Creates the crystal wrapper. """
    from copy import deepcopy
    from .basis import BasisSet
    from .optgeom import OptGeom
    from .electronic import Electronic

    super(Functional, self).__init__()

    self.scf = Electronic()
    """ Holds scf/electronic keywords -- block 3. """
    self.basis   = BasisSet()
    """ Holds definition of basis functions -- block 2. """
    self.optgeom = OptGeom()
    """ Holds definition of geometry optimization -- part of block 1. """
    self.title   = None
    """ Title of the calculation. 
    
        Overriden by the name of the input structure, if it exists.
    """
    self.program = program
    """ Path to crystal program.

        If this attribute is None, then :py:data:`~pylada.crystal_program` is
        used.
    """ 
    self.restart = None
    """ Place holder. """
    if copy is not None:
      for key, value in copy.__dict__.iteritems():
        setattr(self, key, deepcopy(value))

  def __getattr__(self, name):
    """ Pushes scf stuff into instance namespace. """
    from ..error import AttributeError
    if name in self.scf._input: return getattr(self.scf, name)
    raise AttributeError('Unknown attribute {0}.'.format(name))
  def __setattr__(self, name, value):
    """ Pushes scf stuff into instance namespace. """
    from ..tools.input import BaseKeyword
    if isinstance(value, BaseKeyword): 
      if name in ['scf', 'basis', 'optgeom']:
        return super(Functional, self).__setattr__(name, value)
      setattr(self.scf, name, value)
    elif name in self.scf._input:
      setattr(self.scf, name, value)
    else: super(Functional, self).__setattr__(name, value)
  def __delattr__(self, name):
    """ Deletes attributes. """
    if name in self.scf._input: del self.scf._input[name]
    else: super(Functional, self).__delattr__(name)
  def __dir__(self):
    """ List of attributes and members """
    return list( set(self.__dict__.iterkeys()) | set(dir(self.__class__))      \
                 | set(self.scf._input.iterkeys()) )

  def add_keyword(self, name, value=None):
    """ Passes on to :py:attr:`~Functional.scf` """
    return self.scf.add_keyword(name, value)

  def read_input(self, tree, owner=None):
    """ Reads file or string with CRYSTAL input. """
    from ..error import IOError
    from .. import CRYSTAL_geom_blocks as starters
    from .input import find_sym, find_units

    kwargs = {'units': 'angstrom', 'breaksym': True, 'owner': self}

    self.title = tree.keys()[0]
    tree = tree[self.title]
    # read optgeom bit.
    found = False
    for starter in starters:
      if starter in tree.keys(): found = True; break
    if found == False:
      raise IOError('Could not find start of input in file.')
    if 'OPTGEOM' in tree[starter].keys():
      kwargs['breaksym'] = find_sym( tree[starter], 'OPTGEOM', 
                                     kwargs['breaksym'] )
      kwargs['units'] = find_units(tree[starter], 'OPTGEOM', kwargs['units'])
      self.optgeom.read_input(tree[starter]['OPTGEOM'], **kwargs)

    # read basis set
    if 'BASISSET' in tree.keys(): 
      kwargs['breaksym'] = find_sym( tree[starter], 'BASISSET',
                                     kwargs['breaksym'] )
      kwargs['units'] = find_units(tree[starter], 'BASISSET', kwargs['units'])
      self.basis.read_input(tree['BASISSET'], **kwargs)

    # read hamiltonian stuff.
    if 'HAMSCF' in tree.keys():  
      kwargs['breaksym'] = find_sym( tree[starter], 'HAMSCF',
                                     kwargs['breaksym'] )
      kwargs['units'] = find_units(tree[starter], 'HAMSCF', kwargs['units'])
      self.scf.read_input(tree['HAMSCF'], **kwargs)
    if len(self.title) == 0: self.title = None

  def output_map(self, **kwargs):
    """ Dumps CRYSTAL input to string. """
    from ..tools.input import Tree
    from .input import find_units, find_sym

    if 'crystal' not in kwargs: kwargs['crystal'] = self
    root = Tree()
    inner = root

    if 'structure' in kwargs:
      structure = kwargs['structure'].copy()
      # Appends symmetry modification if basis includes a modification of the
      # initial guess.
      modsymm = self.basis.chemod.modisymm(structure)
      if modsymm is not None: structure.append(modsymm)
      # insert name of structure as title.
      if hasattr(structure, 'name'):
        inner = root.descend(structure.name.rstrip().lstrip())
      elif getattr(self, 'title', None) is not None:
        inner = root.descend(self.title.rstrip().lstrip())
      else: inner = root.descend('')
      
      # Output map of the structure
      smap = structure.output_map(**kwargs)
      # To which we add the output map of optgeom
      kwargs['breaksym'] = structure.is_breaksym
      kwargs['units'] = structure.current_units
      smap[0][1].update(self.optgeom.output_map(**kwargs))
      # finally we add that to inner.
      inner.update(smap)
    else: # no structures. Not a meaningful input, but whatever.
      kwargs['breaksym'] = False
      kwargs['units'] = 'angstrom'
      title = getattr(self, 'title', '')
      if title is None: title = ''
      inner = root.descend(title.rstrip().lstrip())
      inner.update(self.optgeom.output_map(**kwargs))

    kwargs['breaksym'] = find_sym(inner, result=kwargs['breaksym'])
    kwargs['units']    = find_units(inner, result=kwargs['units'])
    # now add basis
    inner['BASISSET'] = self.basis.output_map(**kwargs)

    # add scf block
    kwargs['breaksym'] = find_sym(inner, result=kwargs['breaksym'])
    kwargs['units']    = find_units(inner, result=kwargs['units'])
    inner.update(self.scf.output_map(**kwargs))
    # end input and return
    return root

  def print_input(self, **kwargs):
    """ Prints input to string. """
    from .input import print_input
    test = kwargs.pop('test', False)
    map = self.output_map(**kwargs)
    # Does not capitalize input.
    result = map[0][0].rstrip().lstrip() + '\n'
    # Otherwise, everything is standard.
    if test: # If test run, does not actually perform calculation.
      return result + print_input(map[0][1]).rstrip() + '\nTEST\nEND\n'
    return result + print_input(map[0][1]).rstrip() + '\nEND\n'

  def guess_workdir(self, outdir):
    """ Tries and guess working directory. """
    from ..misc import mkdtemp
    from .. import crystal_inplace
    return outdir if crystal_inplace else mkdtemp(prefix='dftcrystal') 

  def bringup(self, structure, outdir, workdir, restart, test):
    """ Creates file environment for run. """
    from os.path import join, abspath, samefile, lexists
    from os import symlink, remove
    from ..misc import copyfile, Changedir
    from ..error import ValueError
    from .. import CRYSTAL_filenames as filenames
    from ..crystal import write
    from .external import External

    # sanity check
    if len(self.basis) == 0:
      raise ValueError('Basis is empty, cannot run CRYSTAL.')
      
    with Changedir(workdir) as cwd:
      # first copies file from current working directory
      if restart is not None: 
        for key, value in filenames.iteritems():
          copyfile( value.format('crystal'), key, nocopyempty=True,
                    symlink=False, nothrow="never" )
      # then copy files from restart.
      if restart is not None:
        for key, value in filenames.iteritems():
          copyfile( join(restart.directory, value.format('crystal')), 
                    key, nocopyempty=True, symlink=False, 
                    nothrow="never" )

      # then creates input file.
      string = self.print_input( crystal=self, structure=structure, 
                                 workdir=workdir, test=test, filework=True )
      with open('crystal.d12', 'w') as file: file.write(string)

    # write additional input to file.
    with Changedir(outdir) as cwd: 
      header = ''.join(['#']*20)
      with open('crystal.out', 'w') as file:
        file.write('{0} {1} {0}\n'.format(header, 'INPUT FILE'))
        file.write(string.rstrip() + '\n')
        file.write('{0} END {1} {0}\n'.format(header, 'INPUT FILE'))
        if isinstance(structure, External):
          file.write('{0} {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
          file.write(write.crystal(structure.initial, None))
          file.write('{0} END {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
        file.write('\n{0} {1} {0}\n'.format(header, 'FUNCTIONAL'))
        file.write(self.__repr__(defaults=False))
        file.write('\n{0} END {1} {0}\n'.format(header, 'FUNCTIONAL'))

    if not samefile(outdir, workdir):
      # Creates symlink to make sure we keep working directory.
      with Changedir(outdir) as cwd:
        with open('crystal.d12', 'w') as file: file.write(string)
        # creates symlink files.
        for filename in ['crystal.err', 'crystal.out', 'ERROR']:
          if lexists(join(workdir, filename)):
            try: remove( join(workdir, filename) )
            except OSError: pass
          try: symlink(abspath(filename), abspath(join(workdir, filename)))
          except OSError: pass
        # adds symlink for unformatted wavefunctions.
        if not lexists(join(workdir, 'fort.9')):
          try: remove(join(workdir, 'fort.9'))
          except OSError: pass
        try: symlink(abspath('crystal.f9'), abspath(join(workdir, 'fort.9')))
        except OSError: pass
        # for optgeom, make sure we bring SCFOUT.LOG 
        if self.optgeom.enabled:
          if self.optgeom.onelog is None or self.optgeom.onelog == False:
            outname = filenames['SCFOUT.LOG'].format('crystal')
            if lexists(join(workdir, 'SCFOUT.LOG')):
              try: remove(join(workdir, 'SCFOUT.LOG'))
              except OSError: pass
            try: 
              symlink( abspath(outname),
                       abspath(join(workdir, 'SCFOUT.LOG')) )
            except OSError: pass
            
        if lexists('workdir'): 
          try: remove('workdir')
          except OSError: pass
        try: symlink(workdir, 'workdir')
        except OSError: pass
    
    # creates a file in the directory, to say we are going to work here
    with open(join(outdir, '.pylada_is_running'), 'w') as file: pass
        

  def bringdown(self, structure, workdir, outdir):
    """ Copies files back to output directory. 
    
        Cats input intO output. Removes workdir if different from outdir
        **and** run was successfull.
    """
    from itertools import chain
    from os import remove
    from os.path import join, samefile, exists
    from shutil import rmtree
    from glob import iglob
    from ..misc import copyfile, Changedir
    from .. import CRYSTAL_filenames, CRYSTAL_delpatterns

    with Changedir(outdir) as cwd:
      # remove 'is running' file marker.
      if exists('.pylada_is_running'):
        try: remove('.pylada_is_running')
        except OSError: pass

      for key, value in CRYSTAL_filenames.iteritems():
        copyfile( join(workdir, key), value.format('crystal'),
                  nocopyempty=True, symlink=False, nothrow="never" )

      header = ''.join(['#']*20)
      if len([0 for filename in iglob(join(workdir, 'ERROR.*'))]):
        string = ""
        for filename in iglob(join(workdir, 'ERROR.*')):
          with open(filename, 'r') as file: string += file.read() + '\n'
        with open('crystal.err', 'w') as out: out.write(string)
        lines = []
        with open('crystal.err', 'r') as out:
          for line in out:
            if len(line.rstrip().lstrip()) == 0: continue
            if line not in lines: lines.append(line.rstrip().lstrip())
        with open('crystal.out', 'a') as out:
          out.write('{0} {1} {0}\n'.format(header, 'ERROR FILE'))
          out.write('\n'.join(lines))
          if len(lines) > 0: out.write('\n')
          out.write('{0} END {1} {0}\n'.format(header, 'ERROR FILE'))
    
    if exists(workdir):
      if samefile(outdir, workdir):
        with Changedir(workdir) as cwd:
          for filepath in chain(*[iglob(u) for u in CRYSTAL_delpatterns]):
            try: remove(filepath)
            except OSError: pass
      elif ExtractBase(outdir).success:
        try: rmtree(workdir)
        except OSError: pass
        try: remove(join(outdir, 'workdir'))
        except OSError: pass
    
    
  
  @stateless
  @assign_attributes(ignore=['overwrite', 'comm', 'workdir', 'test'])
  def iter(self, structure, outdir=None, workdir=None, comm=None,
           overwrite=False, test=False, **kwargs):
    """ Performs a CRYSTAL calculation 
     
        If successfull results (see :py:attr:`extract.Extract.success`) already
        exist in outdir, calculations are not repeated. Instead, an extraction
        object for the stored results are given.

        :param structure:  
            :py:class:`~pylada.crystal.Structure` structure to compute.
        :param outdir:
            Output directory where the results should be stored.  This
            directory will be checked for restart status, eg whether
            calculations already exist. If None, then results are stored in
            current working directory.
        :param comm:
            Holds arguments for executing CRYSTAL externally.
        :param overwrite:
            If True, will overwrite pre-existing results. 
            If False, will check whether a successfull calculation exists. If
            one does, then does not execute. 
        :param test: 
            If True, TEST is inserted at the very end of the run. This makes it
            easier to check that the input is correct.
        :param kwargs:
            Any attribute of the VASP instance can be overidden for
            the duration of this call by passing it as keyword argument.  

        :return: Yields an extractor object if a prior successful run exists.
                 Otherwise, yields a tuple object for executing an external
                 program.

        :note: This functor is stateless as long as self and structure can be
               deepcopied correctly.  

        :raise ExternalRunFailed: when computations do not complete.
        :raise IOError: when outdir exists but is not a directory.
    """ 
    from os import getcwd
    from ..process.program import ProgramProcess
    from ..misc import Changedir, RelativePath
    from .. import crystal_program

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
      self.bringup(structure, outdir, workdir, restart=self.restart, test=test)
      dompi = comm is not None and comm['n'] > 1
      if dompi:
        from ..misc import copyfile
        copyfile('crystal.d12', 'INPUT')

      # figure out the program to launch.
      program = self.program if self.program is not None else crystal_program
      if hasattr(program, '__call__'):
        program = program(self, structure, comm=comm)

      # now creates the process, with a callback when finished.
      onfinish = self.OnFinish(self, structure, workdir, outdir)
      onfail   = self.OnFail(ExtractBase(outdir))
      yield ProgramProcess( program, outdir=workdir, onfinish=onfinish,
                            stdout=None if dompi else ('crystal.out', 'a'), 
                            stderr=('crystal.out', 'a') if dompi else 'crystal.err',
                            stdin=None if dompi else 'crystal.d12', 
                            dompi=dompi, onfail=onfail )
    # yields final extraction object.
    yield ExtractBase(outdir)

  def __call__( self, structure, outdir=None, workdir=None, comm=None,         \
                overwrite=False, test=False, **kwargs):
    from ..error import ExternalRunFailed
    for program in self.iter( structure, outdir=outdir, workdir=workdir,
                              comm=comm, overwrite=overwrite, test=test, 
                              **kwargs ):
      # iterator may yield the result from a prior successfull run. 
      if getattr(program, 'success', False): continue
      # Or may fail return a failed run.
      if not hasattr(program, 'start'): return program
      # otherwise, it should yield a Program tuple to execute.
      program.start(comm)
      program.wait()
    # Last yield should be an extraction object.
    if not program.success:
      raise ExternalRunFailed("CRYSTAL failed to execute correctly.")
    return program
  __call__.__doc__ = iter.__doc__
 
  def __repr__(self, defaults=True, name=None):
    """ Returns representation of this instance """
    from ..tools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    from ..tools.uirepr import template_ui_repr

    results = template_ui_repr(self, imports, name, defaults, ['scf'])
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    scf = self.scf.__ui_repr__(imports, name, getattr(defaults, 'scf', None))
    results.update(scf)
    return results

  def __deepcopy__(self, memo):
    from copy import deepcopy
    result = self.__class__()
    d = deepcopy(self.__dict__, memo)
    result.__dict__.update(d)
    return result

  def __getstate__(self): return self.__dict__
  def __setstate__(self, value):
    self.__dict__.update(value.copy())

  def copy(self): 
    from copy import deepcopy
    return deepcopy(self)

  def test_basis(self, structure, **kwargs):
    """ Returns output of test basis run. """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    try:
      this = deepcopy(self)
      this.basis.add_keyword('TEST')
      tmpdir = mkdtemp()

      result = this(structure, outdir=tmpdir, **kwargs)
      with result.__stdout__() as file: return file.read()
    finally:
      try: rmtree(tmpdir)
      except OSError: pass
  def test_hamiltonian(self, structure, **kwargs):
    """ Returns output of test basis run. """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    try:
      this = deepcopy(self)
      this.add_keyword('TEST')
      tmpdir = mkdtemp()

      result = this(structure, outdir=tmpdir, **kwargs)
      with result.__stdout__() as file: return file.read()
    finally:
      try: rmtree(tmpdir)
      except OSError: pass


  def _ibz(self, structure):
    """ Guesses irreducible k-points from input. """
    from collections import Sequence
    from numpy import abs, floor, dot, any, identity, zeros, array, all
    from numpy.linalg import inv, det
    from ..error import ValueError
    from ..crystal import HFTransform, into_cell
  
    # first figure out supercell and cell in k-space
    cell = structure.eval().cell
    # the supercell is the 
    ksupercell = inv(cell).T
    periodicity = 3
    if abs(cell[2, 2] - 500) < 1e-8: 
      periodicity = 2
      if abs(cell[1, 1] - 500) < 1e-8: 
        periodicity = 1
        if abs(cell[0, 0] - 500) < 1e-8:
          raise ValueError('Non-periodic system.')
  
    mp = self.shrink.mp
    if not isinstance(mp, Sequence):
      mp = [mp] * periodicity
    kcell = ksupercell.copy()
    if periodicity > 2: kcell[:, 2] /= float(mp[2])
    if periodicity > 1: kcell[:, 1] /= float(mp[1])
    if periodicity > 0: kcell[:, 0] /= float(mp[0])
  
    # check those symmetries which leave the reciprocal lattice invariant.
    # also checks for inversion operator.
    symops = []
    invkcell = inv(kcell)
    inverse = -identity(3)
    for op in structure.symmetry_operators:
      if inverse is not None and all(abs(op[:3] - inverse) < 1e-8):
        symops.append(op[:3])
        continue
      transform = dot(invkcell, dot(op[:3], kcell))
      if any(abs(transform - floor(transform)) > 1e-8): continue
      if abs(det(transform) - 1e0) < 1e-8: symops.append(op[:3])
    
    # Now adds inversion operator if it does not exist.
    if inverse is not None:
      for op in [u for u in symops]:
        matrix = dot(op, inverse)
        if all([all(abs(matrix - u) > 1e-8) for u in symops]):
          symops.append(matrix)
        matrix = dot(inverse, op)
        if all([all(abs(matrix - u) > 1e-8) for u in symops]):
          symops.append(matrix)
      if all([all(abs(inverse - u) > 1e-8) for u in symops]):
        symops.append(inverse)
   
    # Now we have a fully consistent symmetry group for the k-mesh.
    # We can actually the work advertised:
    # The Hart-Forcade transforms creates a mapping from any point in a lattice
    # back into a supercell.
    transform = HFTransform(kcell, ksupercell)
    # It can be used to quickly index transformd kpoints and figure out which
    # kpoints are equivalent by symmetry.
    map = zeros(transform.quotient, dtype='bool')
  
    # First we create a list of symmetry operators with the right-hand-side in
    # the index basis and the left-hand-side in cartesian basis.
    invtrans = inv(transform.transform)
    symops = [ dot(op, invtrans) for op in symops 
               if any(abs(op-identity(3)) > 1e-8) ]
  
    # now loop over all k-points an check their symmetries.
    result = []
    for i in xrange(transform.quotient[0]):
      for j in xrange(transform.quotient[1]):
        for k in xrange(transform.quotient[2]):
          i_orig = array([i, j, k])
          if not map[i, j, k]:
            result.append(dot(invtrans, i_orig))
            map[i, j, k] = True
          for op in symops: 
            u, v, w = transform.indices(dot(op, i_orig))
            map[u, v, w] = True
    return into_cell(array(result), ksupercell)

  def _nAOs(self, structure):
    """ Determins number of orbitals. 

        The number of orbitals is *not* reduced by symmetry. 
    """ 
    from ..error import KeyError
    species = [u.type for u in structure.eval()]
    result = 0
    for specie in set(species):
      if specie not in self.basis:
        raise KeyError("Unknown specie {0}.".format(specie))
      dummy = 0
      for shell in self.basis[specie]:
        if shell.type == 's': dummy += 1
        elif shell.type == 'sp': dummy += 4
        elif shell.type == 'p': dummy += 3
        elif shell.type == 'd': dummy += 5
      result += species.count(specie) * dummy
    return result

  def _nb_real_cmplx_kpoints(self, structure):
    """ Number of real and complex k-points. """
    from numpy import dot
    from ..math import is_integer
    # we need to compute real and complex k-points differently.
    # real k-points are those which are on the Brillouin zone edege and Gamma.
    # Their wavefunctions are real through time=reversal and translational
    # symmetry.
    invrecipcell = structure.eval().cell.T
    kpoints = self._ibz(structure)
    kpoints = dot(invrecipcell, kpoints.T).T
    nreal, ncmplx = 0, 0
    for kpoint in kpoints:
      if is_integer(2.*kpoint): nreal += 1
      else: ncmplx += 1
    return nreal, ncmplx

  def mpp_compatible_procs( self, structure, multipleof=1, dof=20,
                            cmplxfac=None ):
    """ Returns list of MPP compatible processor counts.
       
        Processor counts are compatible if:

        - :math:`N_p < N_{AO} * N_k * N_s // dof`
        - :math:`N_p % multipleof == 0`
        - the number of processors per k-point is never prime
        
              Where: 
        
              - :math:`N_p` is the number of cores
              - :math:`N_k` is the number of irreducible kpoints 
        - :math:`N_s` is the number of spins
        - :math:`N_{AO}` is the number of atomic orbitals
        
        The last condition is somewhat more complex and depends upon the total
        number of kpoints, the number of k-points at the Brillouin zone edge,
        and the particular way CRYSTAL_ assigns procs for each k-point. This
        decomposition is computed in :py:meth:`kpoint_blocking`.

        :param structure:
          Crystal structure for which to run the code. 
        :param int multipleof:
	        Pars down results to multiple of this number (eg cores per node).
        :param int dof:
          Minimum degrees of freedom per processor. 
    """
    if cmplxfac is None: cmplxfac = getattr(self, 'cmplxfac', 2) 
    if cmplxfac is None: cmplxfac = 2
    # computes number of real and complex k-points.
    kreal, kcmplx = self._nb_real_cmplx_kpoints(structure)
    # we can now define the max and min number of procs.
    nAOs = self._nAOs(structure)
    # the minimum number of procs should be divisible by multipleof
    weightedprocs = kreal + int(cmplxfac * kcmplx)
    minprocs = int(weightedprocs * 2)
    maxprocs = nAOs * minprocs // dof
    if minprocs % multipleof != 0:
      minprocs = (minprocs // multipleof + 1) * multipleof
 
    # compute all primes, so as to avoid problems later on.
    def allprimes(n):
      """ Computes all primes up to n. """
      primes = list(xrange(2, n+1))
      i = 0
      while i < len(primes):
        fac = primes[i]
        primes = primes[:i+1] + [p for p in primes[i+1:] if p % fac != 0]
        i += 1
      return set([1] + primes)
    primes = allprimes(maxprocs)

    # now loop over potential number of procs and append if ok
    result = []
    for i in xrange(minprocs, maxprocs+1, multipleof):
      blocks = self.kpoint_blocking(structure, i, (kreal, kcmplx))
      if primes.isdisjoint(blocks): result.append(i)
    return result

       
  def kpoint_blocking(self, structure, nprocs, nkpoints=None, cmplxfac=None):
    """ Infers the number of procs per kpoint. 

        This function is helpfull in checking whether MPPcrystal will run
        correcly or not. It returns a list consisting of the number of procs
        per k-point. It is a list since CRYSTAL_ may assign some k-points with
        more processors than others, in order to maximize the number of
        processors used. However, this approach fails when a k-point is
        assigned a prime number of processors, because of some bug somewhere
        (see comment in k_space_MPP.f90: asssing_k_to_proc in CRYSTAL_ code).
        Unfortunately, MPPCrystal does not recover well from this problem.
    """
    # count real and complex k-points in irreducible Brillouin zone.
    if nkpoints is None:
      kreal, kimag = self._nb_real_cmplx_kpoints(structure)
    else: kreal, kimag = nkpoints
    if self.dft.spin: kreal *= 2; kimag *= 2

    # now perform blocking as in crystal.
    # cmplxfac: computational cost of a complex k-point w.r.t. real k-point.
    if cmplxfac is None: cmplxfac = getattr(self, 'cmplxfac', 2)
    if cmplxfac is None: cmplxfac = 2
    # weight: total compuational weight of real+complex k-points.
    weight = cmplxfac*kimag + kreal
    # trivial case where there are more k-points than processors.
    if weight > nprocs: return 0
    # avoids round-off errors
    nreal = int(float(nprocs) / weight)
    nimag = int(cmplxfac * nreal)
    if nreal*kreal + nimag*kimag > nprocs:
      nreal -= 1; nimag = int(cmplxfac * nreal)
      assert nreal*kreal + nimag*kimag <= nprocs
    
    nleft = nprocs - nreal*kreal - nimag*kimag
    if nleft == 0:
      if kreal == 0: return [nimag]
      elif kimag == 0: return [nreal]
      else: return [nreal, nimag]

    result = []
    if kimag != 0:
      imagonly = min(nleft, int(cmplxfac) * kimag)
      nleft -= imagonly
      result.append(nimag + imagonly//kimag)
      if imagonly % kimag != 0: result.append(result[-1] - 1)
      
    if kreal != 0:
      result.append(nreal + nleft // kreal)
      if nleft % kreal != 0: result.append(result[-1] - 1)

    return sorted(set(result))
      
      
    
  def nbelectrons(self, structure):
    """ Number of electrons per formula unit. """
    species = [u.type for u in structure.eval()]
    result = 0
    for specie in set(species):
      if specie not in self.basis:
        raise KeyError("Unknown specie {0}.".format(specie))
      result += sum(u.charge for u in self.basis[species])
    return result
    

  class OnFinish(object):
    """ Called when a run finishes. 
       
        Makes sure files are copied and stuff is pasted at the end of a call to
        CRYSTAL.
    """
    def __init__(self, this, *args):
      self.this = this
      self.args = args
    def __call__(self, *args, **kwargs):
      self.this.bringdown(*self.args)
  class OnFail(object):
    """ Checks whether CRYSTAL run succeeded.

        Crystal reports an error when reaching maximum iteration without
        converging. This screws up how Pylada does things.
    """
    def __init__(self, extract):
      self.extract = extract
    def __call__(self, process, error):
      from ..process import Fail
      self.extract.uncache()
      try: success = self.extract.success
      except: success = False
      if not success:
        raise Fail( 'Crystal failed to run correctly.\n'                       \
                    'It returned with error {0}.'.format(error) )
