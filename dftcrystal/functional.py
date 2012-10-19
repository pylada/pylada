__docformat__ = "restructuredtext en"
from ..tools import stateless, assign_attributes
from .extract import Extract

class Functional(object):
  """ Wrapper for the CRYSTAL program. 
  
      This object provides a pythonic interface to the CRYSTAL_ program. It is
      modelled loosely after CRYSTAL_'s input:

        - The OPTGEOM keyword of the first code block can be accessed through
          the :py:class:`optgeom <lada.dftcrystal.optgeom.OptGeom>`
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
          <lada.dftcrystal.basis.BasisSet>` :py:attr:`attribute <basis>`. It
          allows to set the basis set itself, as well as keywords specific to
          the basis set:

          .. code-block:: python
            
            functional.basis['H'] = [Shell('s', a0=(15.0, 1.0), a1=(8.0, 1.0)), ...]
            functional.basis.ghosts = [1, 4, 6], False

        - The third input block (Hamiltonian and miscellaneous) can be accessed
          directly through the functional, or, alternatively, *via* the
          :py:class:`scf <lada.dftcrystal.electronic.Electronic>`
          :py:attr:`attribute <scf>`.
        
          .. code-block:: python
        
            functional.scf.dft.b3lyp = True
            functional.dft.b3lyp = True
        
            functional.tolinteg = [8] * 4 + [14]
        
          The first two lines are exactly equivalent. The third line could also
          be written as ``functional.scf.tolinteg = ...``.

  """
  Extract = Extract
  """ Extraction class. """
  __ui_name__ = 'functional'
  """ Name used in user-friendly representation """
  def __init__(self, copy=None, program=None, **kwargs):
    """ Creates the crystal wrapper. """
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

        If this attribute is None, then :py:data:`~lada.crystal_program` is
        used.
    """ 
    self.restart = None
    """ Place holder. """

  def __getattr__(self, name):
    """ Pushes scf stuff into instance namespace. """
    from ..error import ValueError
    if name in self.scf._input: return getattr(self.scf, name)
    raise ValueError('Unknown attribute {0}.'.format(name))
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

    self.title = tree.keys()[0]
    tree = tree[self.title]
    # read optgeom bit.
    found = False
    for starter in starters:
      if starter in tree.keys(): found = True; break
    if found == False:
      raise IOError('Could not find start of input in file.')
    if 'OPTGEOM' in tree[starter].keys():
      self.optgeom.read_input(tree[starter]['OPTGEOM'], owner=self)

    # read basis set
    if 'BASISSET' in tree.keys(): 
      self.basis.read_input(tree['BASISSET'], owner=self)

    # read hamiltonian stuff.
    if 'HAMSCF' in tree.keys():  
      self.scf.read_input(tree['HAMSCF'], owner=self)

  def output_map(self, **kwargs):
    """ Dumps CRYSTAL input to string. """
    from ..tools.input import Tree

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
      smap[0][1].update(self.optgeom.output_map(**kwargs))
      # finally we add that to inner.
      inner.update(smap)
    else: # no structures. Not a meaningful input, but whatever.
      title = getattr(self, 'title', '')
      if title is None: title = ''
      inner = root.descend(title.rstrip().lstrip())
      inner.update(self.optgeom.output_map(**kwargs))

    # now add basis
    inner['BASISSET'] = self.basis.output_map(**kwargs)

    # add scf block
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
    from os import environ
    from tempfile import mkdtemp
    from datetime import datetime
    from .. import crystal_inplace
    if crystal_inplace: return outdir
    rootdir = environ.get('PBS_TMPDIR', outdir)
    return mkdtemp( prefix='dftcrystal_{0!s}'.format(datetime.today()), 
                    dir=rootdir )

  def bringup(self, structure, outdir, workdir, restart, test):
    """ Creates file environment for run. """
    from os.path import join, abspath, samefile, lexists
    from os import symlink, remove
    from ..misc import copyfile, Changedir
    from ..error import ValueError
    from .. import CRYSTAL_filenames as filenames

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

    with Changedir(outdir) as cwd: pass
    if not samefile(outdir, workdir):
      # Creates symlink to make sure we keep working directory.
      with Changedir(outdir) as cwd:
        with open('crystal.d12', 'w') as file: file.write(string)
        with open('crystal.out', 'w') as file: pass
        with open('crystal.err', 'w') as file: pass
        with open('ERROR', 'w') as file: pass
        # creates symlink files.
        for filename in ['crystal.err', 'crystal.out', 'ERROR']:
          if lexists(join(workdir, filename)):
            try: remove( join(workdir, filename) )
            except: pass
          symlink(abspath(filename), abspath(join(workdir, filename)))
        # for optgeom, make sure we bring SCFOUT.LOG 
        if self.optgeom.enabled:
          if self.optgeom.onelog is None or self.optgeom.onelog == False:
            outname = filenames['SCFOUT.LOG'].format('crystal')
            with open(outname, 'w') as file: pass
            if lexists(join(workdir, 'SCFOUT.LOG')):
              try: remove(join(workdir, 'SCFOUT.LOG'))
              except: pass
            try: 
              symlink( abspath(outname),
                       abspath(join(workdir, 'SCFOUT.LOG')) )
            except: pass
            
        if lexists('workdir'): 
          try: remove('workdir')
          except: pass
        try: symlink(workdir, 'workdir')
        except: pass
    
    # creates a file in the directory, to say we are going to work here
    with open(join(outdir, '.lada_is_running'), 'w') as file: pass
        

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
    from .external import External
    from ..crystal import write
    from ..misc import copyfile, Changedir
    from .. import CRYSTAL_filenames, CRYSTAL_delpatterns

    with Changedir(outdir) as cwd:
      for key, value in CRYSTAL_filenames.iteritems():
        copyfile( join(workdir, key), value.format('crystal'),
                  nocopyempty=True, symlink=False, nothrow="never" )

      with open('crystal.d12', 'r') as file: input = file.read()
      with open('crystal.out', 'r') as file: output = file.read()
      header = ''.join(['#']*20)
      with open('crystal.out', 'w') as file:
        file.write('{0} {1} {0}\n'.format(header, 'INPUT FILE'))
        input = input.rstrip()
        if input[-1] != '\n': input += '\n'
        file.write(input)
        file.write('{0} END {1} {0}\n'.format(header, 'INPUT FILE'))
        if isinstance(structure, External):
          file.write('{0} {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
          file.write(write.crystal(structure.initial, None))
          file.write('{0} END {1} {0}\n'.format(header, 'INITIAL STRUCTURE'))
        file.write(output)
        file.write('\n{0} {1} {0}\n'.format(header, 'FUNCTIONAL'))
        file.write(self.__repr__(defaults=False))
        file.write('\n{0} END {1} {0}\n'.format(header, 'FUNCTIONAL'))
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
          out.write('{0} END {1} {0}\n'.format(header, 'ERROR FILE'))

      # remove 'is running' file marker.
      if exists('.lada_is_running'):
        try: remove('.lada_is_running')
        except: pass
    
    if samefile(outdir, workdir):
      with Changedir(workdir) as cwd:
        for filepath in chain(*[iglob(u) for u in CRYSTAL_delpatterns]):
          try: remove(filepath)
          except: pass
    elif Extract(outdir).success:
      try: rmtree(workdir)
      except: pass
      try: remove(join(outdir, 'workdir'))
      except: pass
    
    
  
  @assign_attributes(ignore=['overwrite', 'comm', 'workdir'])
  @stateless
  def iter(self, structure, outdir=None, workdir=None, comm=None,
           overwrite=False, test=False, **kwargs):
    """ Performs a CRYSTAL calculation 
     
        If successfull results (see :py:attr:`extract.Extract.success`) already
        exist in outdir, calculations are not repeated. Instead, an extraction
        object for the stored results are given.

        :param structure:  
            :py:class:`~lada.crystal.Structure` structure to compute.
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

        :raise RuntimeError: when computations do not complete.
        :raise IOError: when outdir exists but is not a directory.
    """ 
    from os.path import abspath
    from os import getcwd
    from ..process.program import ProgramProcess
    from ..misc import Changedir
    from .. import crystal_program

    # check for pre-existing and successfull run.
    if not overwrite:
      extract = self.Extract(outdir)
      if extract.success:
        yield extract # in which case, returns extraction object.
        return
    
    if outdir == None: outdir = getcwd()
    if workdir == None: workdir = self.guess_workdir(outdir)

    outdir = abspath(outdir)
    workdir = abspath(workdir)
    with Changedir(workdir) as tmpdir: 

      # writes/copies files before launching.
      self.bringup(structure, outdir, workdir, restart=self.restart, test=test)
      dompi = comm is not None
      if dompi:
        from ..misc import copyfile
        copyfile('crystal.d12', 'INPUT')

      # figure out the program to launch.
      program = self.program if self.program is not None else crystal_program
      if hasattr(program, '__call__'):
        program = program(self, structure, comm=comm)

      # now creates the process, with a callback when finished.
      onfinish = self.OnFinish(self, structure, workdir, outdir)
      onfail   = self.OnFail(outdir)
      yield ProgramProcess( program, outdir=workdir, onfinish=onfinish,
                            stdout=None if dompi else 'crystal.out', 
                            stderr='crystal.out' if dompi else 'crystal.err',
                            stdin=None if dompi else 'crystal.d12', 
                            dompi=dompi, onfail=onfail )
    # yields final extraction object.
    yield Extract(outdir)

  def __call__( self, structure, outdir=None, workdir=None, comm=None,         \
                overwrite=False, test=False, **kwargs):
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
      raise RuntimeError("CRYSTAL failed to execute correctly.")
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
      except: pass
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
      except: pass


  def _ibz(self, structure, **kwargs):
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

  def _nAOs(self, structure, **kwargs):
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

  def ndogs(self, structure, **kwargs):
    """ Number of degrees of liberty. 

        This is the main number of MPPcrystal's golden rule for choosing the
        number of cores: :math:`N_{cores} = \\frac{N_{spin}N_{AO}N_{IBZ}}{20}`.
        This routine returns the numerator in the fraction.
    """
    return (2 if self.dft.spin else 1)                                         \
           * self._nAOs(structure, **kwargs)                                   \
           * len(self._ibz(structure, **kwargs))

  def nbelectrons(self, structure, **kwargs):
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
        converging. This is screws up hwo LaDa does things.
    """
    def __init__(self, outdir):
      self.outdir = outdir
    def __call__(self, process, error):
      from ..process import Fail
      try: success = Extract(self.outdir).success
      except: success = False
      if not success:
        raise Fail( 'Crystal failed to run correctly.\n'                       \
                    'It returned with error {0}.'.format(error) )
