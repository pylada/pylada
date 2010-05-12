""" Interface module for pescan. """
from sys import exit
from ..opt import __load_pescan_in_global_namespace__
if __load_pescan_in_global_namespace__:
  from DLFCN import RTLD_NOW as _RTLD_NOW, RTLD_GLOBAL as _RTLD_GLOBAL
  from sys import getdlopenflags as _getdlopenflags, setdlopenflags as _setdlopenflags
  flags = _getdlopenflags()
  _setdlopenflags(_RTLD_NOW|_RTLD_GLOBAL)
  import _escan
  _setdlopenflags(flags)
else: import _escan
from ..opt.decorators import add_setter, broadcast_result
from bandstructure import band_structure
import _extract
Extract = _extract.Extract

def _is_in_sync(comm, which = [0]):
  from boost.mpi import broadcast
  if comm == None: return 
  print "sync ", comm.rank, which[0]
  which[0] += 1
  m = broadcast(comm, "666" if comm.rank == 0 else None, 0)
  return m == "666"


def call_escan(comm, atom="atom_input", pot="pot_input", escan="escan_input"):
  """ Calls escan functional in current directory.

      Before calling the functional, the files are propagated such that each
      proc can read its own. What is more, the inputs are modified such that
      L{pot} does refer to L{atom}, and L{escan} to both atomic config and
      potential files.
      @param comm: Processes on which to execute.
      @type comm: boost.mpi.Communicator or mpi4py
      @param atom: file with atomic input (from vff). 
      @param pot: file with input to potential generation.
      @param escan: file with input to escan itself.
  """
  from os import remove
  from boost.mpi import world, broadcast
  from _escan import _call_escan, _call_genpot
  
  # escan is free for all, with every proc reading at the same time.
  # hence it must be from different files.
  atominput = "%s.%i" % (atom, world.rank)
  potinput = "%s.%i" % ("pot.input", world.rank)
  potoutput = "%s.%i" % ("pot.output", world.rank)
  escaninput = "%s.%i" % ("escan_input", world.rank)

  # propagates atomic input first
  buffer = None
  if comm.rank == 0:
    with open(atom, "r") as file: buffer = file.read()
  buffer = broadcast(comm, buffer, 0)
  with open(atominput, "w") as file: file.write(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(pot, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[0], atom )
  buffer = broadcast(comm, buffer, 0)
  with open(potinput, "w") as file: file.writelines(buffer)

  # propagates pot input second
  buffer = None
  if comm.rank == 0:
    with open(escan, "r") as file: buffer = file.readlines()
    buffer[0] = buffer[0].replace( buffer[0].split()[1], potoutput )
    buffer[12] = buffer[12].replace( buffer[12].split()[1], atominput )
  buffer = broadcast(comm, buffer, 0)
  with open(escaninput, "w") as file: file.writelines(buffer)

  #  calls escan at last. 
  _call_genpot(comm)
  _call_escan(comm)

  remove(atominput)
  remove(potinput)
  remove(escaninput)

folded_spectrum = 0
""" Folded spectrum method. """
all_electron = 1
""" All electron method. """
localH = 0
""" Local hamiltonian. """
nonlocalH = 1
""" Non-local hamiltonian. """
soH = 2
""" Spin-orbit hamiltonian. """

class AtomicPotential(object):
  """ Holds parameters to atomic potentials. """
  def __init__(self, filepath, nonlocal=None, s=None, p=None, d=None, pnl=None, dnl=None):
    from os.path import abspath

    self.filepath = abspath(filepath)
    """ Path to pseudopotential file. """
    self.nonlocal = None if nonlocal == None else abspath(nonlocal)
    """ Path to non-local part, or None. """
    self.s =  s if s != None else 0
    """ s parameter """
    self.p =  p if p != None else 0
    """ p parameter """
    self.d =  d if d != None else 0
    """ d parameter """
    self.pnl = pnl if pnl != None else 0
    """ p non-local parameter """
    self.dnl = dnl if dnl != None else 0
    """ d non-local parameter """

  def __repr__(self):
    """ Tuple representation of self. """
    result = "\"%s\"" % (self.filepath)
    if self.nonlocal == None: result += ", None"
    else: result += ", \"%s\"" % (self.nonlocal)
    result += ", %f, %f, %f, %f, %f" % (self.s, self.p, self.d, self.pnl, self.dnl)
    return result

  @broadcast_result(key=True)
  def get_izz(self, comm = None):
    """ Returns izz string greped from pseudopotential file. """
    with open(self.filepath, "r") as file:
      return file.readline().split()[1]


class Escan(object):
  """ Performs PESCAN calculations, from structure relaxation to wavefunctions. """

  def __init__(self, workdir=None):
    """ Initializes a PESCAN functional. """
    from numpy import zeros
    from ..vff import Vff

    super(Escan, self).__init__()
    self.vff = Vff() 
    """ The L{Vff} functional with which to relax a structure. """
    self.OUTCAR = "escan_out" 
    """ Escan output file. """
    self.ERRCAR = "escan_err"
    """ Escan error file. """
    self.WAVECAR = "wavefunctions"
    """ Wavefunction file (in g-space). """
    self.eref = None
    """ Reference energy for folded spectrum method.
    
        Set to None for all electron diagonalization.
    """
    self.cutoff = 8.2
    """ Cutoff energy for plane-wave expansion. """
    self.smooth = 1
    """ Smooth potential scaling. """
    self.kinetic_scaling = 1
    """ Smooth kinetic energy scaling. """
    self.nbstates = 1
    """ Number of states to compute. """
    self.itermax = 20
    """ Maximum number of iterations when minimzing wavefuntions. """
    self.nlines = 50
    """ Conjugate gradient line minimization parameter. """
    self.tolerance = 1e-12
    """ Requested minimization accuracy. """
    self.input_wavefunctions = []
    """ Indices of wavefunctions to read in. """
    self.INWAVECAR = "input_escan_wavefunctions"
    """ Filename of input wavefunctions. """
    self.kpoint = zeros((3,1), dtype="float64")
    """ k-point at which to perform calculations. """
    self.potential = soH
    """ Type of hamiltonian to use. """
    self.rspace_cutoff = 5
    """ real-space projector cutoff. """
    self.atomic_potentials = None
    """ Parameters to atomic potentials. """
    self.fft_mesh = (18, 18, 18)
    """ Fourrier Transform mesh. """
    self.dnc_mesh = None
    """ Divide and conquer mesh.
     
        Set to None if no divide and conquer is required. 
    """
    self.overlap_mesh = None
    """ Overlap of divide and conquer mesh. 
     
        Set to None if no divide and conquer is required. 
    """
    self.do_relax = True
    """ If true, the structure is relaxed using vff.
          
        Otherwise, the file L{_POSCAR} from the output directory is used. 
    """
    self.do_genpot = True
    """ If true, the potential is generated.
          
        Otherwise, the file L{_POTCAR} from the output directory is used. 
    """
    self.do_escan = True
    """ If true, calculations are performed. """

    self._POSCAR = "atomic_input"
    """ Private reference to the atomic input file. """
    self._workdir = workdir
    """ Private reference to the working directory. """
    self._POTCAR = "pot.output"
    """ Private reference to the potential file generated by genpot. """
    self._maskr = "maskr"
    """ Private reference to the maskr projector file. """
    self._INCAR = "escan_input"
    """ Private reference to the escan input file. """
    self._GENCAR = "pot.input"
    """ Private reference to the genpot input file. """


  def _get_workdir(self): return self._workdir
  def _set_workdir(self, workdir):
    from os.path import abspath
    self._workdir = abspath(workdir) if workdir != None else None
  workdir = property( _get_workdir, _set_workdir, 
                      """ Working directory where calculations are performed. 
                      
                          Absolute path is determined when workdir is set.
                      """ )

  def _get_maskr(self): return self._maskr
  def _set_maskr(self, workdir): 
    from os.path import abspath
    self._maskr = abspath(workdir) if workdir != None else None
  maskr = property( _get_maskr, _set_maskr, 
                      """ maskr file for projectors.
                      
                          Absolute path is determined when maskr is set.
                      """
                    )
  def _get_lattice(self): return self.vff.lattice
  def _set_lattice(self, arg): self.vff.lattice = arg
  lattice = property( _get_lattice, _set_lattice, "Lattice to use with escan and vff.")

  def _add_potential(self, args):
    """ Adds an atomic potential to escan. """
    assert len(args) > 2, RuntimeError("Atomic  potentials need at least two parameters.")
    assert len(args) < 9, RuntimeError("Too many parameters when setting atomic potentials.")
    if self.atomic_potentials == None: self.atomic_potentials = []
    self.atomic_potentials.append( AtomicPotential(*args) )
  add_potential = add_setter( _add_potential,\
                              """ Adds atomic potential to escan.
                                  
                                  - first argument is the path to the atomic potential. 
                                    The absolute path is deduced when set.
                                  - second argument is the path to the non-local potential file. 
                                    If None, then no non-local argument is added. Defaults to None.
                                  - third trough seventh arguments are the
                                    non-local potential parameters s, p, d,
                                    pnl, dnl. Defaults to None (eg 0).
                              """ )

  def __repr__(self):
    result  = str(self.vff)
    result += "# Escan definition.\n"
    result += "escan = %s()\n" % (self.__class__.__name__)
    result += "escan.eref                  = %s\n"\
              % ( "None" if self.eref == None else str(self.eref) )
    result += "escan.cutoff                = %f\n" % (self.cutoff)
    result += "escan.smooth                = %f\n" % (self.smooth)
    result += "escan.kinetic_scaling       = %f\n" % (self.kinetic_scaling)
    result += "escan.nbstates              = %i\n" % (self.nbstates)
    result += "escan.itermax               = %i\n" % (self.itermax)
    result += "escan.nlines                = %i\n" % (self.nlines)
    result += "escan.tolerance             = %e\n" % (self.tolerance)
    result += "escan.rspace_cutoff         = %f\n" % (self.rspace_cutoff)
    result += "escan.fft_mesh              = %i, %i, %i\n" % self.fft_mesh
    result += "escan.do_genpot             = %s\n" % ("True" if self.do_genpot else "False")
    result += "escan.do_escan              = %s\n" % ("True" if self.do_escan else "False")
    result += "escan.do_relax              = %s\n" % ("True" if self.do_relax else "False")
    result += "escan.input_wavefunctions   = %s\n" % (repr(self.input_wavefunctions))
    result += "escan.kpoint                = %s\n" % (repr(self.kpoint))
    result += "escan.dnc_mesh              = %s\n" % (repr(self.dnc_mesh))
    result += "escan.overlap_mesh          = %s\n" % (repr(self.overlap_mesh))
    if self.potential == localH:
      result += "escan.potential             = localH\n"
    elif self.potential == nonlocalH:
      result += "escan.potential             = localH\n"
    elif self.potential == soH:
      result += "escan.potential             = soH\n"
    else: raise RuntimeError("unknown hamiltonnian %i." % (soH))
    for pot in self.atomic_potentials:
      result += "escan.add_potential         = %s\n" % (repr(pot))
    result += "escan.INWAVECAR = \"%s\"\n" % (self.INWAVECAR)
    result += "escan.ERRCAR = \"%s\"\n" % (self.ERRCAR)
    result += "escan.WAVECAR = \"%s\"\n" % (self.WAVECAR)
    result += "escan.workdir = \"%s\"\n" % (self.workdir)
    result += "escan.maskr = \"%s\"\n" % (self.maskr)
    result += "escan._INCAR = \"%s\"\n" % (self._INCAR)
    result += "escan._POTCAR = \"%s\"\n" % (self._POTCAR)
    result += "escan._GENCAR = \"%s\"\n" % (self._GENCAR)
    result += "# End of escan definition."
    return result

  def __call__(self, structure, outdir = None, comm = None, **kwargs):
    """ Performs calculation """
    from copy import deepcopy
    from os import getcwd
    from os.path import exists, isdir, abspath, basename, join
    from shutil import copyfile
    from boost.mpi import world, broadcast
    from ..opt.changedir import Changedir

    if comm == None: comm = world
    if outdir == None: outdir = getcwd()

    # make this functor stateless.
    this      = deepcopy(self)
    outdir    = abspath(outdir)

    # if other keyword arguments are present, then they are assumed to be
    # attributes of self, with value to use for calculations launch. 
    for key in kwargs.keys():
      if hasattr(this, key): setattr(this, key, kwargs[key])
      else: raise NameError( "Escan does not have an %s attribute." % (key) )

    # checks if outdir contains a successful run.
    if broadcast(comm, exists(outdir) if comm.rank == 0 else None, 0):
      extract = Extract(comm = comm, directory = outdir, escan = this)
      if extract.success: return extract # in which case, returns extraction object.
      comm.barrier() # makes sure directory is not created by other proc!

    this._run(structure, outdir, comm)

    with Changedir(outdir) as cwd:
      for file in  [ this._POSCAR + "." + str(comm.rank),\
                     this._cout(comm) if this._cout(comm) != "/dev/null" else None,\
                     this._cerr(comm) if this._cerr(comm) != "/dev/null" else None,\
                     this.vff._cout(comm) if this.vff._cout(comm) != "/dev/null" else None,\
                     this.vff._cerr(comm) if this.vff._cerr(comm) != "/dev/null" else None,\
                     this.WAVECAR if comm.rank == 0  else None ]:
        if file == None: continue
        copyfile( join(this.workdir, basename(file)), basename(file) )

    return Extract(comm = comm, directory = outdir, escan = this)

  def _cout(self, comm):
    """ Creates output name. """
    if self.OUTCAR == None: return "/dev/null"
    return self.OUTCAR if comm.rank == 0 else self.OUTCAR + "." + str(comm.rank)

  def _cerr(self, comm):
    """ Creates error name. """
    if self.OUTCAR == None: return "/dev/null"
    return self.ERRCAR if comm.rank == 0 else self.ERRCAR + "." + str(comm.rank)


  def _run(self, structure, outdir, comm):
    """ Performs escan calculation. """
    import time
    from os.path import join
    from ..opt.changedir import Changedir

    if self.do_genpot == False and self.do_relax == False and self.do_escan == False:
      print "Nothing to do? no relaxation, no genpot, no escan?" 
      return None
    timing = time.time() 
    local_time = time.localtime() 

    # prints some output first
    cout, cerr = self._cout(comm), self._cerr(comm)
    with Changedir(self.workdir) as cwd:
      with open(cout, "w") as file: 
        print >>file, "# Escan calculation on ", time.strftime("%m/%d/%y", local_time),\
                      " at ", time.strftime("%I:%M:%S %p", local_time)
        if len(structure.name) != 0: print "# Structure named ", structure.name 
        print >>file, repr(self)
        print >>file, "# Performing calculations. "
      
      # makes calls to run
      self._run_vff(structure, "", comm, cout)
      if not (self.do_genpot and self.do_escan): return None 
      self._run_genpot(comm, cout, cerr)
      if not self.do_escan: return None
      self._run_escan(comm, structure.scale, cout, cerr)
     
      with open(cout, "a") as file: 
        timing = time.time() - timing
        hour = int(float(timing/3600e0))
        minute = int(float((timing - hour*3600)/60e0))
        second = (timing - hour*3600-minute*60)
        extract = Extract(comm = comm, directory = outdir, escan = self)
        assert not extract.success, RuntimeError("Escan calculations did not complete.")
        print >> file, "# Computed ESCAN in: %i:%i:%f."  % (hour, minute, second) 

  def _run_vff(self, structure, outdir, comm, cout):
    """ Gets atomic input ready, with or without relaxation. """
    from shutil import copyfile
    from os.path import join, samefile

    poscar = self._POSCAR + "." + str(comm.rank)
    if not self.do_relax:
      if exists(join(outdir, poscar)):
        copyfile(join(outdir, poscar), join(self.workdir, poscar))
      vff.write_escan_input(poscar, structure)
      return 

    out = self.vff(structure, outdir=outdir, comm=comm)
    assert out.success, RuntimeError("VFF relaxation did not succeed.")
    out.write_escan_input(poscar, out.structure)

    # copies vff output to stdout. This way, only one outcar.
    if comm.rank == 0 and samefile(out.OUTCAR, cout):
      with open(out.OUTCAR) as file_in: 
        with open(cout) as file_out: 
          for line in file_in:
            if line.find("# VFF calculation on ") != -1: print >>file_out, line[:-1]
            if line == "# Performing VFF calculations. ": break
          print >>file_out, line[:-1]
          for line in file_in:
            if line.find("# Computed VFF in:") != -1: break
            print >>file_out, line[:-1]
          print >>file_out, line[:-1]


  def _run_genpot(self, comm, cout, cerr):
    """ Runs genpot only """
    from ._escan import _call_genpot
    from shutil import copyfile
    from os.path import basename
    from ..opt import redirect

    if not self.do_genpot: 
      if exists(join(outdir, self._GENCAR)):
        copyfile(join(outdir, self._GENCAR), join(self.workdir, self._GENCAR))
        return
      else:
        with open(cout, "a") as file:
          print >> file, "Could not find genpot potential %s."\
                         "Recomputing." % (join(outdir, self._GENCAR))

    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    dnc_mesh = self.dnc_mesh if self.dnc_mesh != None else self.fft_mesh
    overlap_mesh = self.overlap_mesh if self.overlap_mesh != None else (0,0,0)
    with open(self._GENCAR + "." + str(comm.rank), "w") as file:
      file.write( "%s\n%i %i %i\n%i %i %i\n%i %i %i\n%f\n%i\n"\
                  % ( self._POSCAR, self.fft_mesh[0], self.fft_mesh[1], self.fft_mesh[2], \
                      dnc_mesh[0], dnc_mesh[1], dnc_mesh[2],\
                      overlap_mesh[0], overlap_mesh[1], overlap_mesh[2], self.cutoff,\
                      len(self.atomic_potentials) ))
      for pot in self.atomic_potentials:
        # adds to list of potentials
        filepath = basename(pot.filepath)
        file.write(filepath + "\n") 
        # copy potential files as well.
        if comm.rank == 0: copyfile(pot.filepath, filepath)
        if pot.nonlocal != None and comm.rank == 0: 
          copyfile(pot.nonlocal, basename(pot.nonlocal))
    if comm.rank == 0: copyfile(self.maskr, basename(self.maskr))

    comm.barrier() # syncs all procs
    with redirect(fout=cout, ferr=cerr, append=True) as oestreams: 
      _call_genpot(comm)

  def _run_escan(self, comm, scale, cout, cerr):
    """ Runs escan only """
    from shutil import copyfile
    from os.path import basename
    from numpy.linalg import norm
    from ._escan import _call_escan
    from ..opt import redirect
    from ..physics import a0

    assert self.atomic_potentials != None, RuntimeError("Atomic potentials are not set.")
    # Creates temporary input file and creates functional
    with open(self._INCAR + "." + str(comm.rank), "w") as file:
      print >> file, "1 %s.%i" % (self._POTCAR, comm.rank) 
      print >> file, "2 %s" % (self.WAVECAR) 
      print >> file, "3 %i # %s" % ((1, "folded spectrum") if self.eref != None else (2, "all electron"))
      print >> file, "4 %f %f %f %f # Eref, cutoff, smooth, kinetic scaling"\
                     % ( self.eref if self.eref != None else 0,\
                         self.cutoff, self.smooth, self.kinetic_scaling )
      if self.potential != soH or norm(self.kpoint) < 1e-6: 
        print >> file, "5 %i # number of states" % (self.nbstates/2)
      else: print >> file, "5 %i # number of states" % (self.nbstates)

      print >> file, "6 %i %i %e # itermax, nllines, tolerance" % (self.itermax, self.nlines, self.tolerance)
      nowfns = self.input_wavefunctions == None
      if not nowfns: nowfns = len(self.input_wavefunctions) == 0
      if nowfns: print >> file, "7 0 # no input wfns\n8 0 # wfns indices"
      else:
        print >> file, "7 %i\n8 %i" % (len(self.input_wavefunctions), self.input_wavefunctions[0])
        for u in self.input_wavefunctions[1:]:
          print >> file, str(u),
        print >> file
      print >> file, "9 %s # input wavefunction filename" % (self.INWAVECAR)

      print >> file, "10 0 1 1 1 0"

      if norm(self.kpoint) < 1e-12: print >> file, "11 0 0 0 0 0"
      else: print >> file, "11 1 %i %i %i" % tuple((self.kpoint * scale / a0("A")).flat)
      
      if   self.potential == localH: print >> file, "12 1 # local hamiltonian" 
      elif self.potential == nonlocalH: print >> file, "12 2 # non-local hamiltonian" 
      elif self.potential == soH: print >> file, "12 3 # spin orbit hamiltonian" 
      else: raise RuntimeError("Unknown potential requested.")
      
      print >> file, "13 ", self._POSCAR + "." + str(comm.rank)
      print >> file, "14 ", self.rspace_cutoff, "# real-space cutoff" 

      if self.potential != soH:
        print >> file, "15 ", 0, "# Number of spin-orbit potentials"
      else:
        print >> file, "15 ", len(self.atomic_potentials), "# Number of spin-orbit potentials"
        for i, pot in enumerate(self.atomic_potentials):
          filepath = basename(pot.nonlocal)
          print >> file, i + 16, filepath, pot.get_izz(comm), pot.s , pot.p, pot.d, pot.pnl, pot.dnl

    comm.barrier() # syncs all procs to make sure we are reading from same file.
    if comm.rank == 0: copyfile(self.maskr, basename(self.maskr))
    with redirect(fout=cout, ferr=cerr, append=True) as oestreams: 
      _call_escan(comm)
