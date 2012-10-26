""" Python wrapper for the CRYSTAL dft code. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Extract', 'AffineTransform', 'DisplaceAtoms', 'InsertAtoms',
            'Marker', 'ModifySymmetry', 'RemoveAtoms', 'Slabcut', 'read',
            'Slabinfo','Crystal', 'Molecule', 'Slab', 'Functional', 'Shell',
            'read_gaussian_basisset', 'MassExtract', 'read_input', 'exec_input',
            'Elastic', 'External', 'Supercell', 'Supercon' ]

from .basis import Shell
from .functional import Functional
from .extract import Extract, MassExtract
from .geometry import AffineTransform, DisplaceAtoms, InsertAtoms,              \
                      ModifySymmetry, RemoveAtoms, Slabcut, Slabinfo, Slab,     \
                      Elastic, Supercell, Supercon
from .crystal import Crystal
from .external import External
from .molecule import Molecule
registered = { 'atomrot':    AffineTransform,
               'atomdisp':   DisplaceAtoms,
               'atominse':   InsertAtoms,
               'modisymm':   ModifySymmetry,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'elastic':    Elastic,
               'supercel':   Supercell,
               'supercon':   Supercon,
               'crystal':    Crystal,
               'external':   External,
               'molecule':   Molecule }
""" Keywords used in creating the geometry. """

def read(path):
  """ Reads CRYSTAL input from file

      Reads a file and finds the CRYSTAL_ input within. It then parses the file
      and creates a :py:class:`~functional.Functional` and
      :py:class:`~molecule.Molecule` (or derived) object.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      :returns: A two-tuple containing the structure and the functional

      :rtype: (:py:class:`~molecule.Molecule`, :py:class:`~functional.Functional`)

      .. _CRYSTAL: http://www.crystal.unito.it/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from .parse import parse
  b = parse(path)
  key = b.keys()[0]

  func = Functional()
  crys = Crystal()
  func.read_input(b)
  crys.read_input(b[key]['CRYSTAL'])
  crys.name = key
  return crys, func


def read_gaussian_basisset(path):
  """ Reads a GAUSSIAN94 basis set

      This is meant to read input files from the `basis set exchange website`__
      in the GAUSSIAN94_ format.

      :param path:
        Can be one of the following:

          - a file object obtained from open_
          - a string without '\\n' in which case it should be a path
          - a string with '\\n', in whcih case it should be the input itself
          - an iterator over the lines of an input file

      .. __: https://bse.pnl.gov/bse/portal
      .. _GAUSSIAN94: http://www.gaussian.com/
      .. _open: http://docs.python.org/library/functions.html#open
  """
  from ..error import GrepError
  from ..misc import RelativePath
  if isinstance(path, str): 
    if path.find('\n') == -1:
      with open(RelativePath(path).path) as file: return read_gaussian_basisset(file)
    else:
      return read_gaussian_basisset(path.split('\n').__iter__())

  for line in path: 
    if set(['*']) == set(line.rstrip().lstrip()): break

  # read specie type
  try: line = path.next().split()
  except StopIteration: raise GrepError('Unexpected end of file')
  specie = line[0]
  result = {specie: []}
  
  try: 
    while True:
      line = path.next().split()
      if len(line) != 3:
        line = path.next().split()
        if len(line) != 2: break
        specie = line[0]
        result[specie] = []
        continue
      type, n, scale = line[0].lower(), int(line[1]), float(line[2])
      shell = Shell(type)
      for i in xrange(n): shell.append(*path.next().split())
      result[specie].append(shell)

  except StopIteration: pass

  return result

def read_input(filepath="input.py", namespace=None):
  """ Specialized read_input function for dftcrystal. 
  
      :Parameters: 
        filepath : str
          A path to the input file.
        namespace : dict
          Additional names to include in the local namespace when evaluating
          the input file.

      It add a few names to the input-file's namespace. 
  """
  from ..misc import read_input

  # names we need to create input.
  input_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input':
      input_dict[k] = globals()[k]
  if namespace is not None: input_dict.update(namespace)
  return read_input(filepath, input_dict)

def exec_input( script, global_dict=None, local_dict=None,
                paths=None, name=None ):
  """ Specialized exec_input function for vasp. """
  from ..misc import exec_input

  # names we need to create input.
  if global_dict is None: global_dict = {}
  for k in __all__:
    if k != 'read_input' and k != 'exec_input':
      global_dict[k] = globals()[k]
  return exec_input(script, global_dict, local_dict, paths, name)

from lada import is_interactive
if is_interactive:
  from IPython.core.magic import register_line_magic

  @register_line_magic
  def complete_crystal(line):
    """ Adds input file to beginning of output file where necessary.

        Works only for dftcrystal and should only be enabled when the
        dftcrystal is loaded from ipython.

        A complete file should contain all the information necessary to recreate
        that file. Unfortunately, this is generally not the case with CRYSTAL's
        standard output, at least not without thorough double-guessing.

        This function adds the .d12 file if it not already there, as well as the
        input structure if it is in "external" format.
    """
    from lada import is_interactive
    if not is_interactive: return
    ip = get_ipython() # should be defined in interactive mode.
    if 'collect' not in ip.user_ns:
      print 'No jobfolder loaded in memory.'
      return
    collect = ip.user_ns['collect']
    jobparams = ip.user_ns['jobparams']
    for name, extractor in collect.iteritems():
      if not hasattr(extractor, '_complete_output'): continue

      structure = jobparams[name].structure
      if extractor._complete_output(structure):
        print 'Corrected ', extractor.directory
  del complete_crystal

  @register_line_magic
  def crystal_mppcount(line):
    """ Number of processors on which group of jobs can be sent.

	It can be a bit difficult to determine whether an MPPcrystal job will
	run correctly (e.g. with k-point parallelization), especially if more
	than one job is involved and each job must have the same number of
        processors.
     
	This magic function determines likely core counts (per job) across a
	job-folder. It can take a single integer argument to filter the
	processor counts to multiples of that integer (e.g. to run on complete
        nodes only). 
    """
    # figure multipleof
    shell = get_ipython()
    line = line.rstrip().lstrip()
    if len(line) == 0: multipleof = 1
    elif line in shell.user_ns: multipleof = shell.user_ns[line]
    else: multipleof = eval(line, shell.user_ns.copy())

    # figure jobparams
    if 'jobparams' not in shell.user_ns:
      print 'No jobfolder loaded in memory.'
      return
    jobparams = shell.user_ns['jobparams']

    allsets = []
    for job in jobparams['/'].itervalues():
      if job.is_tagged: continue
      dummy = job.functional.mpp_compatible_procs( job.structure,
                                                   multipleof=multipleof )
      allsets.append(dummy)
    
    maxprocs = min(max(m) for m in allsets)
    result = set(list(xrange(maxprocs)))
    for s in allsets: result &= set(s)
    return sorted(result)

  del crystal_mppcount
  del register_line_magic

del is_interactive

