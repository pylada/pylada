""" Subpackage containing extraction methods for VASP. """
__docformat__  = 'restructuredtext en'
__all__ = ['ExtractBase']
from quantities import g, cm, eV
from ...functools import make_cached
from ...functools.extract import search_factory

OutcarSearchMixin = search_factory('OutcarSearchMixin', 'OUTCAR', __name__)

class ExtractBase(object):
  """ Implementation class for extracting data from VASP output """

  def __init__(self):
    """ Initializes the extraction class. """
    super(ExtractBase, self).__init__()

  @property 
  @make_cached
  def ialgo(self):
    """ Returns the kind of algorithms. """
    result = self._find_first_OUTCAR(r"""^\s*IALGO\s*=\s*(\d+)\s*""")
    return int(result.group(1))

  @property 
  @make_cached
  def algo(self):
    """ Returns the kind of algorithms. """
    return { 68: 'Fast', 38: 'Normal', 48: 'Very Fast', 58: 'Conjugate',
             53: 'Damped', 4: 'Subrot', 90: 'Exact', 2: 'Nothing'}[self.ialgo]

  @property
  def is_dft(self):
    """ True if this is a DFT calculation, as opposed to GW. """
    try: return self.algo not in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
    except: return False
  @property
  def is_gw(self):
    """ True if this is a GW calculation, as opposed to DFT. """
    try: return self.algo in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
    except: return False
    
  @property 
  def encut(self):
    """ Energy cutoff. """
    return float(self._find_first_OUTCAR(r"ENCUT\s*=\s*(\S+)").group(1)) * eV


  @property
  @make_cached
  def functional(self):
    """ Returns vasp functional used for calculation.

        Requires the functional to be pasted at the end of the calculations. 
    """
    from re import compile, M
    from ...misc import exec_input
    regex = compile('#+ FUNCTIONAL #+\n((.|\n)*)\n#+ END FUNCTIONAL #+')
    with self.__outcar__() as file: result = regex.search(file.read())
    if result is None: return None
    return exec_input(result.group(1)).functional

  @property
  @make_cached
  def success(self):
    """ Checks that VASP run has completed. 

        At this point, checks for the existence of OUTCAR.
        Then checks that timing stuff is present at end of OUTCAR.
    """
    from os.path import exists, join

    if hasattr(self, "OUTCAR"): 
      for path in [self.OUTCAR]:
        if not exists(join(self.directory, path)): return False
      
    regex = r"""General\s+timing\s+and\s+accounting\s+informations\s+for\s+this\s+job"""
    return self._find_last_OUTCAR(regex) is not None

  @property
  @make_cached
  def datetime(self):
    """ Greps execution date and time. """
    from datetime import datetime
    regex = r"""executed on .*\n"""
    result = self._find_first_OUTCAR(regex)
    if result is None: return
    result = result.group(0)
    result = result[result.find('date')+4:].rstrip().lstrip()
    return datetime.strptime(result, '%Y.%m.%d  %H:%M:%S')


  @property
  def initial_structure(self):
    """ Structure at start of calculations. """
    from re import compile, M
    from numpy import array, zeros, dot
    from numpy.linalg import inv
    from ...crystal import Structure
    from ...misc import exec_input

    try:
      regex = compile('#+ INITIAL STRUCTURE #+\n((.|\n)*)\n#+ END INITIAL STRUCTURE #+')
      with self.__outcar__() as file: result = regex.search(file.read(), M)
      if result is not None: return exec_input(result.group(1)).structure
    except: pass

    result = Structure()
    with self.__outcar__() as file: 
      atom_index, cell_index = None, None
      cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
      atom_re = compile(r"""^\s*position\s+of\s+ions\s+in\s+fractional\s+coordinates""")
      for line in file:
        if cell_re.search(line) is not None: break
      data = []
      for i in range(3):
        data.append(file.next().split())
      try: 
        for i in range(3): result.cell[:,i] = array(data[i][:3], dtype='float64')
      except: 
        for i in range(3): result.cell[i, :] = array(data[i][-3:], dtype='float64')
        cell = inv(cell)
      for line in file:
        if atom_re.search(line) is not None: break
      for specie, n in zip(self.species,self.stoichiometry):
        for i, line in zip(range(n), file):
          data = line.split()
          result.add_atom(pos=dot(result.cell, array(data, dtype='float64')), type=specie)
    return result

  @property 
  def _catted_contcar(self):
    """ Returns contcar found at end of OUTCAR. """
    from re import compile
    from StringIO import StringIO
    from ...crystal import read
    with self.__outcar__() as file: lines = file.readlines()
    begin_contcar_re = compile(r"""#+\s+CONTCAR\s+#+""")
    end_contcar_re = compile(r"""#+\s+END\s+CONTCAR\s+#+""")
    start, end = None, None
    for i, line in enumerate(lines[::-1]):
      if begin_contcar_re.match(line) is not None: start = -i; break;
      if end_contcar_re.match(line) is not None: end = -i
    if start is None or end is None:
      raise IOError("Could not find catted contcar.")
    stringio = StringIO("".join(lines[start:end if end != 0 else -1]))
    return read.poscar(stringio, self.species)

  @property
  @make_cached
  def _grep_structure(self):
    """ Greps cell and positions from OUTCAR. """
    from re import compile
    from ...crystal import Structure
    from numpy.linalg import inv
    from numpy import array

    with self.__outcar__() as file: lines = file.readlines()
    atom_index, cell_index = None, None
    atom_re = compile(r"""^\s*POSITION\s+""")
    cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
    for index, line in enumerate(lines[::-1]):
      if atom_re.search(line) is not None: atom_index = index - 1
      if cell_re.search(line) is not None: cell_index = index; break
    assert atom_index is not None and cell_index is not None,\
           RuntimeError("Could not find structure description in OUTCAR.")
    result = Structure()
    try: 
      for i in range(3):
        result.cell[:,i] = array(lines[-cell_index+i].split()[:3], dtype="float64")
    except: 
      for i in range(3): result.cell[i,:] = array(lines[-cell_index+i].split()[-3:], dtype="float64")
      result.cell = inv(result.cell)
    species = [type for type, n in zip(self.species, self.stoichiometry) for i in xrange(n)]
    while atom_index > 0 and len(lines[-atom_index].split()) == 6:
      result.add_atom( pos=array(lines[-atom_index].split()[:3], dtype="float64"), 
                       type=species.pop(-1) )
      atom_index -= 1

    return result

  @property
  @make_cached
  def structure(self):
    """ Greps structure and total energy from OUTCAR. """
    from quantities import eV

    if self.nsw == 0 or self.ibrion == -1:
      return self.starting_structure

    try: result = self._catted_contcar;
    except:
      result = self._contcar_structure
      try: result = self._contcar_structure
      except:
        result = self._grep_structure

    # tries to find adequate name for structure.
    try: name = self.system
    except RuntimeError: name = ''
    if len(name) == 0 or name == 'POSCAR created by SUPERPOSCAR':
      try: title = self.system
      except RuntimeError: title = ''
      if len(title) != 0: result.name = title
    else: result.name = name
    
    if self.is_dft: result.energy = self.total_energy 

    try: initial = self.initial_structure
    except: pass
    else:
      for key, value in initial.__dict__.iteritems():
        if not hasattr(result, key): setattr(result, key, value)
      for a, b in zip(result, initial):
        for key, value in b.__dict__.iteritems():
          if not hasattr(a, key): setattr(a, key, value)
    # adds magnetization.
    try: magnetization = self.magnetization
    except: pass
    else:
      if magnetization is not None:
        for atom, mag in zip(result, magnetization): atom.magmom = sum(mag)
    # adds stress.
    try: result.stress = self.stress
    except: pass
    # adds forces.
    try: forces = self.forces
    except: pass
    else:
      for force, atom in zip(forces, result): atom.force = force
    return result

  @property
  @make_cached
  def LDAUType(self):
    """ Type of LDA+U performed. """
    type = self._find_first_OUTCAR(r"""LDAUTYPE\s*=\s*(\d+)""")
    if type == None: return None
    type = int(type.group(1))
    if type == 1: return "liechtenstein"
    elif type == 2: return "dudarev"
    return type

  @property
  @make_cached
  def HubbardU_NLEP(self):
    """ Hubbard U/NLEP parameters. """
    from ..specie import U as ldaU, nlep
    from re import M
    type = self._find_first_OUTCAR(r"""LDAUTYPE\s*=\s*(\d+)""")
    if type == None: return {}
    type = int(type.group(1))

    species = tuple([ u.group(1) for u in self._search_OUTCAR(r"""VRHFIN\s*=\s*(\S+)\s*:""") ])

    # first look for standard VASP parameters.
    groups = r"""\s*((?:(?:\+|-)?\d+(?:\.\d+)?\s*)+)\s*\n"""
    regex = r"""\s*angular\s+momentum\s+for\s+each\s+species\s+LDAUL\s+=""" + groups \
          + r"""\s*U\s+\(eV\)\s+for\s+each\s+species\s+LDAUU\s+="""         + groups \
          + r"""\s*J\s+\(eV\)\s+for\s+each\s+species\s+LDAUJ\s+="""         + groups
    result = {} 
    for k in species: result[k] = []
    for found in self._search_OUTCAR(regex, M):
      moment = found.group(1).split()
      LDAU   = found.group(2).split()
      LDAJ   = found.group(3).split()
      for specie, m, U, J in zip(species, moment, LDAU, LDAJ):
        if int(m) != -1: result[specie].append(ldaU(type, int(m), float(U), float(J)))
    for key in result.keys():
      if len(result[key]) == 0: del result[key]
    if len(result) != 0: return result

    # then look for nlep parameters.
    regex = r"""\s*angular\s+momentum\s+for\s+each\s+species,\s+LDAU#\s*(?:\d+)\s*:\s*L\s*=""" + groups \
          + r"""\s*U\s+\(eV\)\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*U\s*=""" + groups \
          + r"""\s*J\s+\(eV\)\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*J\s*=""" + groups \
          + r"""\s*nlep\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*O\s*=""" + groups 
    result = {} 
    for k in species: result[k] = []
    for found in self._search_OUTCAR(regex, M):
      moment = found.group(1).split()
      LDAU   = found.group(2).split()
      LDAJ   = found.group(3).split()
      funcs  = found.group(4).split()
      for specie, m, U, J, func in zip(species, moment, LDAU, LDAJ, funcs):
        if int(m) == -1: continue
        if int(func) == 1: result[specie].append(ldaU(type, int(m), float(U), float(J)))
        else: result[specie].append(nlep(type, int(m), float(U), float(J)))
    for key in result.keys():
      if len(result[key]) == 0: del result[key]
    return result

  @property
  @make_cached
  def pseudopotential(self):
    """ Title of the first POTCAR. """
    return self._find_first_OUTCAR(r"""POTCAR:.*""").group(0).split()[1]


  @property
  @make_cached
  def volume(self): 
    """ Unit-cell volume. """
    from numpy import abs
    from numpy.linalg import det
    from quantities import angstrom
    return abs(det(self.structure.scale * self.structure.cell)) * angstrom**3

  @property 
  @make_cached
  def reciprocal_volume(self):
    """ Reciprocal space volume (including 2pi). """
    from numpy import abs, pi
    from numpy.linalg import det, inv
    from quantities import angstrom
    return abs(det(inv(self.structure.scale * self.structure.cell))) * (2e0*pi/angstrom)**3

  @property
  @make_cached
  def density(self):
    """ Computes density of the material. """
    from quantities import atomic_mass_unit as amu
    from ... import periodic_table as pt
    result = 0e0 * amu;
    for atom, n in zip(self.species, self.stoichiometry): 
      if atom not in pt.__dict__: return None;
      result += pt.__dict__[atom].atomic_weight * float(n) * amu 
    return (result / self.volume).rescale(g/cm**3)

  @property
  @make_cached
  def _contcar_structure(self):
    """ Greps structure from CONTCAR. """
    from ...crystal import read
    from quantities import eV

    result = read.poscar(self.__contcar__(), self.species)
    result.energy = float(self.total_energy.rescale(eV)) if self.is_dft else 0e0
    return result

  @property
  @make_cached
  def stoichiometry(self):
    """ Stoichiometry of the compound. """
    result = self._find_first_OUTCAR(r"""\s*ions\s+per\s+type\s*=.*$""")
    if result is None: return None
    return [int(u) for u in result.group(0).split()[4:]]
  @property
  def ions_per_specie(self): 
    """ Alias for stoichiometry. """
    return self.stoichiometry

  @property
  @make_cached
  def species(self):
    """ Greps species from OUTCAR. """
    return [ u.group(1) for u in self._search_OUTCAR(r"""VRHFIN\s*=\s*(\S+)\s*:""") ]

  @property
  @make_cached
  def isif(self):
    """ Greps ISIF from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISIF\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))
  
  @property
  @make_cached
  def nsw(self):
    """ Greps NSW from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*NSW\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def ismear(self):
    """ Greps smearing function from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISMEAR\s*=\s*(\d+);""")
    if  result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def sigma(self):
    """ Greps smearing function from OUTCAR. """
    from numpy import array
    from quantities import eV
    result = self._find_first_OUTCAR(r"""\s*ISMEAR\s*=\s*(?:\d+)\s*;\s*SIGMA\s*=\s+(.*)\s+br""")
    if  result is None: return None
    result = result.group(1).rstrip().lstrip().split()
    if len(result) == 1: return float(result[0]) * eV
    return array(result, dtype="float64") * eV
  
  @property
  def relaxation(self): 
    """ Returns the kind of relaxation performed in this calculation. """
    from ..incar import Relaxation
    result = Relaxation(self.isif)
    if self.nsw != 50: result.nsw = self.nsw
    if self.ibrion != 2 or abs(self.potim - 0.5) > 1e-8: result.ibrion = self.ibrion
    if abs(self.potim - 0.5) > 1e-8: result.potim = self.potim
    return result.value
  
  @property
  def smearing(self): 
    """ Returns the kind of smearing performed in this calculation. """
    from quantities import eV
    from ..incar import Smearing
    result = Smearing()
    result.ismear = self.ismear
    result.sigma = self.sigma if abs(self.sigma - 0.2*eV) > 1e-8 else None
    return result.value

  @property
  @make_cached
  def ibrion(self):
    """ Greps IBRION from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*IBRION\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def potim(self):
    """ Greps POTIM from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*POTIM\s*=\s*(-?\S+)\s+""")
    if result is None: return None
    return float(result.group(1))

  @property
  @make_cached
  def lorbit(self):
    """ Greps LORBIT from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*LORBIT\s*=\s*(\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def isym(self):
    """ Greps ISYM from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISYM\s*=\s*(\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def nupdown(self):
    """ Greps NUPDOWN from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*NUPDOWN\s*=\s*(\S+)\s+""")
    if result is None: return None
    return float(result.group(1))

  @property
  @make_cached
  def lmaxmix(self):
    """ Greps LMAXMIX from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*LMAXMIX\s*=\s*(\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def istart(self):
    """ Greps ISTART from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISTART\s*=\s*(\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def icharg(self):
    """ Greps ICHARG from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ICHARG\s*=\s*(\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @make_cached
  def precision(self):
    """ Greps PREC from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*PREC\s*=\s*(\S*)\s+""")
    if result is None: return None
    if result.group(1) == "accura": return "accurate"
    return result.group(1)

  @property
  @make_cached
  def ediff(self):
    """ Greps EDIFF from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*EDIFF\s*=\s*(\S+)\s+""")
    if result is None: return None
    return float(result.group(1))

  @property
  @make_cached
  def ediffg(self):
    """ Greps EDIFFG from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*EDIFFG\s*=\s*(\S+)\s+""")
    if result is None: return None
    return float(result.group(1))

  @property
  @make_cached
  def kpoints(self):
    """ Greps k-points from OUTCAR.
    
        Numpy array where each row is a k-vector in cartesian units. 
    """
    from re import compile
    from numpy import array

    result = []
    found_generated = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
    found_read = compile(r"""k-points in units of 2pi/SCALE and weight""")
    with self.__outcar__() as file:
      found = 0
      for line in file:
        if found_generated.search(line) is not None: found = 1; break
        elif found_read.search(line) is not None: found = 2; break
      if found == 1:
        found = compile(r"""Following\s+cartesian\s+coordinates:""")
        for line in file:
          if found.search(line) is not None: break
        file.next()
        for line in file:
          data = line.split()
          if len(data) != 4: break;
          result.append( data[:3] )
        return array(result, dtype="float64") 
      if found == 2:
        for line in file:
          data = line.split()
          if len(data) == 0: break
          result.append(data[:3])
        return array(result, dtype="float64") 

  @property
  @make_cached
  def multiplicity(self):
    """ Greps multiplicity of each k-point from OUTCAR. """
    from re import compile
    from numpy import array

    result = []
    # case where kpoints where generated by vasp.
    found_generated = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
    # case where kpoints where given by hand.
    found_read = compile(r"""k-points in units of 2pi/SCALE and weight""")
    with self.__outcar__() as file:
      found = 0
      for line in file:
        if found_generated.search(line) is not None: found = 1; break
        elif found_read.search(line) is not None: found = 2; break
      if found == 1:
        found = compile(r"""Following\s+cartesian\s+coordinates:""")
        for line in file:
          if found.search(line) is not None: break
        file.next()
        for line in file:
          data = line.split()
          if len(data) != 4: break;
          result.append( data[-1] )
        return array(result, dtype="float64") 
      elif found == 2:
        for line in file:
          data = line.split()
          if len(data) == 0: break
          result.append(float(data[3]))
        return array([round(r*float(len(result))) for r in result], dtype="float64")

  @property 
  @make_cached
  def ispin(self):
    """ Greps ISPIN from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*ISPIN\s*=\s*(1|2)\s+""")
    assert result is not None, RuntimeError("Could not extract ISPIN from OUTCAR.")
    return int(result.group(1))

  @property
  @make_cached
  def name(self):
    """ Greps POSCAR title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*POSCAR\s*=.*$""")
    assert result is not None, RuntimeError("Could not extract POSCAR title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:]
    return result.rstrip().lstrip()

  @property
  @make_cached
  def system(self):
    """ Greps system title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*SYSTEM\s*=.*$""")
    assert result is not None, RuntimeError("Could not extract SYSTEM title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:].rstrip().lstrip()
    if result[0] == '"': result = result[1:]
    if result[-1] == '"': result = result[:-1]
    return result

  def _unpolarized_values(self, which):
    """ Returns spin-unpolarized eigenvalues and occupations. """
    from re import compile, finditer
    import re

    with self.__outcar__() as file: lines = file.readlines()
    # Finds last first kpoint.
    spin_comp1_re = compile(r"\s*k-point\s+1\s*:\s*(\S+)\s+(\S+)\s+(\S+)\s*")
    found = None
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found is not None: break
    assert found is not None, RuntimeError("Could not extract eigenvalues/occupation from OUTCAR.")

    # now greps actual results.
    if self.is_dft:
      kp_re = r"\s*k-point\s+(?:(?:\d|\*)+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+band\s+energies\s+occupation\s*\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 2, 3
    else: 
      kp_re = r"\s*k-point\s+(?:(?:\d|\*)+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+.*\n\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)"\
              r"\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 3, 8
    results = []
    for kp in finditer(kp_re, "".join(lines[-i-1:]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results.append([float(u[which]) for u in dummy if len(u) == cols])
    return results

  def _spin_polarized_values(self, which):
    """ Returns spin-polarized eigenvalues and occupations. """
    from re import compile, finditer
    import re

    with self.__outcar__() as file: lines = file.readlines()
    # Finds last spin components.
    spin_comp1_re = compile(r"""\s*spin\s+component\s+(1|2)\s*$""")
    spins = [None,None]
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found is None: continue
      if found.group(1) == '1': 
        assert spins[1] is not None, \
               RuntimeError("Could not find two spin components in OUTCAR.")
        spins[0] = i
        break
      else:  spins[1] = i
    assert spins[0] is not None and spins[1] is not None,\
           RuntimeError("Could not extract eigenvalues/occupation from OUTCAR.")

    # now greps actual results.
    if self.is_dft:
      kp_re = r"\s*k-point\s+(?:(?:\d|\*)+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+band\s+energies\s+occupation\s*\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 2, 3
    else: 
      kp_re = r"\s*k-point\s+(?:(?:\d|\*)+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+.*\n\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)"\
              r"\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 3, 8
    results = [ [], [] ]
    for kp in finditer(kp_re, "".join(lines[-spins[0]:-spins[1]]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results[0].append([float(u[which]) for u in dummy if len(u) == cols])
    for kp in finditer(kp_re, "".join(lines[-spins[1]:]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results[1].append([u[which] for u in dummy if len(u) == cols])
    return results

  @property
  @make_cached
  def ionic_charges(self):
    """ Greps ionic_charges from OUTCAR."""
    regex = """^\s*ZVAL\s*=\s*(.*)$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find ionic_charges in OUTCAR")
    return [float(u) for u in result.group(1).split()]

  @property
  @make_cached
  def valence(self):
    """ Greps number of valence bands from OUTCAR."""
    ionic = self.ionic_charges
    species = self.species
    atoms = [u.type for u in self.structure]
    result = 0
    for c, s in zip(ionic, species): result += c * atoms.count(s)
    return result
  
  @property
  @make_cached
  def nelect(self):
    """ Greps nelect from OUTCAR."""
    regex = """^\s*NELECT\s*=\s*(\S+)\s+total\s+number\s+of\s+electrons\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) 

  @property
  def extraelectron(self):
    """ Returns charge state of the system. """
    return self.nelect - self.valence

  @property
  def nonscf(self):
    """ True if non-selfconsistent calculation. """
    return self.icharg >= 10

  @property
  def lwave(self):
    """ Greps LWAVE from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*LWAVE\s*=\s*(\S)""")
    if result == None: return None
    return result.group(1) == 'T' 
   
  @property
  def lcharg(self):
    """ Greps LWAVE from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*LCHARG\s*=\s*(\S)""")
    if result == None: return None
    return result.group(1) == 'T' 

  @property
  def lvtot(self):
    """ Greps LVTOT from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*LVTOT\s*=\s*(\S)""")
    if result == None: return None
    return result.group(1) == 'T' 

  @property
  def nelm(self):
    """ Greps NELM from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*NELM\s*=\s*(\d+)""")
    if result == None: return None
    return int(result.group(1))

  @property
  def nelmdl(self):
    """ Greps NELMDL from OUTCAR. """
    regex = r"""^\s*NELM\s*=\s*\d+\s*;"""\
             """\s*NELMIN\s*=\s*\d+\s*;"""\
             """\s*NELMDL\s*=\s*(-?\d+)"""
    result = self._find_first_OUTCAR(regex)
    if result == None: return None
    return int(result.group(1))

  @property
  def nelmin(self):
    """ Greps NELMIN from OUTCAR. """
    regex = r"""^\s*NELM\s*=\s*\d+\s*;\s*NELMIN\s*=\s*(\d+)"""
    result = self._find_first_OUTCAR(regex)
    if result == None: return None
    return int(result.group(1))

  @property
  @make_cached
  def nbands(self):
    """ Number of bands in calculation. """
    result = self._find_first_OUTCAR("""NBANDS\s*=\s*(\d+)""")
    assert result is not None, RuntimeError("Could not find NBANDS in OUTCAR.")
    return int(result.group(1))

  @property
  @make_cached
  def nbprocs(self):
    """ Number of bands in calculation. """
    result = self._find_first_OUTCAR("""running\s+on\s+(\d+)\s+nodes""")
    assert result is not None, RuntimeError("Could not find number of processes in OUTCAR.")
    return int(result.group(1))


  def iterfiles(self, **kwargs):
    """ iterates over input/output files. 
    
        :kwarg errors: Include stderr files.
        :type errors: bool
        :kwarg incar: Include INCAR file
        :type incar: bool
        :kwarg wavecar: Include WAVECAR file
        :type wavecar: bool
        :kwarg doscar: Include CHGCAR file
        :type doscar: bool
        :kwarg chgcar: Include CHGCAR file
        :type chgcar: bool
        :kwarg poscar: Include POSCAR file
        :type poscar: bool
        :kwarg contcar: Include CONTCAR file
        :type contcar: bool
        :kwarg procar: Include PROCAR file
        :type procar: bool
    """
    from os.path import exists, join
    from glob import iglob
    from itertools import chain
    files = [self.OUTCAR]
    try: files.append(self.functional.STDOUT)
    except: pass
    if kwargs.get('errors', False): 
      try: files.append(self.functional.STDERR)
      except: pass
    if kwargs.get('incar', False):   files.append('INCAR')
    if kwargs.get('wavecar', False): files.append('WAVECAR')
    if kwargs.get('doscar', False):  files.append('DOSCAR')
    if kwargs.get('chgcar', False):  files.append('CHGCAR')
    if kwargs.get('poscar', False):  files.append('POSCAR')
    if kwargs.get('contcar', False): files.append('CONTCAR')
    if kwargs.get('procar', False): files.append('PROCAR')
    for file in files:
      file = join(self.directory, file)
      if exists(file): yield file
    # Add RelaxCellShape directories.
    for dir in chain( iglob(join(self.directory, "relax_cellshape/[0-9]/")),
                      iglob(join(self.directory, "relax_cellshape/[0-9][0-9]/")),
                      iglob(join(self.directory, "relax_ions/[0-9]/")),
                      iglob(join(self.directory, "relax_ions/[0-9][0-9]/")) ):
      a = self.__class__(dir)
      for file in a.iterfiles(**kwargs): yield file

  @property
  @make_cached
  def energy_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find sigma0 energy in OUTCAR")
    return float(result.group(2)) * eV

  @property
  @make_cached
  def energies_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(2)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @property
  @make_cached
  def all_total_energies(self):
    """ Greps total energies for all electronic steps from OUTCAR."""
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy =\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @property
  def fermi0K(self):
    """ Fermi energy at zero kelvin.
    
        This procedure recomputes the fermi energy using a step-function. 
        It avoids negative occupation numbers. As such, it may be different
        from the fermi energy given by vasp, depeneding on the smearing and the
        smearing function.
    """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from operator import itemgetter
    from numpy import rollaxis, sum

    if self.ispin == 2: eigs = rollaxis(self.eigenvalues, 0, 2)
    else: eigs = self.eigenvalues
    array = [(e, m) for u, m in zip(eigs, self.multiplicity) for e in u.flat]
    array = sorted(array, key=itemgetter(0))
    occ = 0e0
    factor = (2.0 if self.ispin == 1 else 1.0)/float(sum(self.multiplicity))
    for i, (e, w) in enumerate(array):
      occ += w*factor
      if occ >= self.valence - 1e-8: return e * self.eigenvalues.units
    return None

  @property
  def halfmetallic(self):
    """ True if the material is half-metallic. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import max, abs
    if self.ispin == 1: return False

    if self.cbm - self.vbm > 0.01 * self.cbm.units: return False

    fermi = self.fermi0K
    vbm0 = max(self.eigenvalues[0][self.eigenvalues[0]<=float(fermi)+1e-8]) 
    vbm1 = max(self.eigenvalues[1][self.eigenvalues[1]<=float(fermi)+1e-8])
    return abs(float(vbm0-vbm1)) > 2e0 * min(float(fermi - vbm0), float(fermi - vbm1))

  @property
  def cbm(self):
    """ Returns Condunction Band Minimum. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import min
    return min(self.eigenvalues[self.eigenvalues>self.fermi0K+1e-8*self.fermi0K.units])

  @property
  def vbm(self):
    """ Returns Valence Band Maximum. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import max
    from quantities import eV
    return max(self.eigenvalues[self.eigenvalues<=self.fermi0K+1e-8*self.fermi0K.units])

  @property
  @make_cached
  def total_energies(self):
    """ Greps total energies for all ionic steps from OUTCAR."""
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @property
  @make_cached
  def total_energy(self):
    """ Greps total energy from OUTCAR."""
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import eV
    regex = """energy\s+without\s+entropy=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) * eV

  energy = total_energy
  """ Alias for total_energy. """

  @property
  @make_cached
  def fermi_energy(self):
    """ Greps fermi energy from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import eV
    regex = r"""E-fermi\s*:\s*(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find fermi energy in OUTCAR")
    return float(result.group(1)) * eV

  @property
  @make_cached
  def moment(self):
    """ Returns magnetic moment from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find magnetic moment in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
  def pressures(self):
    """ Greps all pressures from OUTCAR """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find pressures in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find pressures in OUTCAR")
    return result * kB

  @property
  @make_cached
  def pressure(self):
    """ Greps last pressure from OUTCAR """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find pressure in OUTCAR")
    return float(result.group(1)) * kB

  @property
  @make_cached
  def alphabet(self):
    """ Greps alpha+bet from OUTCAR """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find alpha+bet in OUTCAR")
    return float(result.group(3))

  @property
  @make_cached
  def xc_g0(self):
    """ Greps XC(G=0) from OUTCAR """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find xc(G=0) in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
  def pulay_pressure(self):
    """ Greps pressure from OUTCAR """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find pulay pressure in OUTCAR")
    return float(result.group(2)) * kB

  @property
  @make_cached
  def fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from re import search, M as re_M
    from numpy import array
    regex = r"""(?!I would recommend the setting:\s*\n)"""\
             """\s*dimension x,y,z NGX =\s+(\d+)\s+NGY =\s+(\d+)\s+NGZ =\s+(\d+)"""
    with self.__outcar__() as file: result = search(regex, file.read(), re_M)
    assert result is not None, RuntimeError("Could not FFT grid in OUTCAR.""")
    return array([int(result.group(1)), int(result.group(2)), int(result.group(3))])

  @property
  @make_cached
  def recommended_fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from re import search, M as re_M
    from numpy import array
    regex = r"""I would recommend the setting:\s*\n"""\
             """\s*dimension x,y,z NGX =\s+(\d+)\s+NGY =\s+(\d+)\s+NGZ =\s+(\d+)"""
    with self.__outcar__() as file: result = search(regex, file.read(), re_M)
    assert result is not None, RuntimeError("Could not recommended FFT grid in OUTCAR.""")
    return array([int(result.group(1)), int(result.group(2)), int(result.group(3))])

  def _get_partial_charges_magnetization(self, grep):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included. Implementation also used for magnetization.
    """
    import re 
    from numpy import array

    result = []
    with self.__outcar__() as file: lines = file.readlines()
    found = re.compile(grep) 
    for index in xrange(1, len(lines)+1):
      if found.search(lines[-index]) is not None: break 
    if index == len(lines): return None
    index -= 4
    line_re = re.compile(r"""^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$""")
    for i in xrange(0, index): 
      match = line_re.match(lines[-index+i])
      if match is None: break
      result.append([float(match.group(j)) for j in range(1, 5)])
    return array(result, dtype="float64")

  @property
  @make_cached
  def partial_charges(self):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    return self._get_partial_charges_magnetization(r"""\s*total\s+charge\s*$""")

  @property
  @make_cached
  def magnetization(self):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    return self._get_partial_charges_magnetization(r"""^\s*magnetization\s*\(x\)\s*$""")

  @property
  @make_cached
  def eigenvalues(self):
    """ Greps eigenvalues from OUTCAR. 

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    if not (self.is_dft or self.is_gw):
      raise AttributeError('Neither a DFT not a GW calculation.')
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(1), dtype="float64") * eV

  @property
  @make_cached
  def occupations(self):
    """ Greps occupations from OUTCAR. 

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    if not (self.is_dft or self.is_gw):
      raise AttributeError('Neither a DFT not a GW calculation.')
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(2), dtype="float64")
    return array(self._unpolarized_values(2), dtype="float64") 

  @property
  @make_cached
  def qp_eigenvalues(self):
    """ Greps quasi-particle eigenvalues from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    if not self.is_gw: raise AttributeError('not a GW calculation.')
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(2), dtype="float64") * eV

  @property
  @make_cached
  def self_energies(self):
    """ Greps self-energies of each eigenvalue from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    if not self.is_gw: raise AttributeError('not a GW calculation.')
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(3), dtype="float64") * eV

  @property
  @make_cached
  def electropot(self):
    """ Greps average atomic electrostatic potentials from OUTCAR. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from re import compile, X as reX
    from numpy import array
    from quantities import eV

    with self.__outcar__() as file: lines = file.readlines()
    regex = compile(r"""average\s+\(electrostatic\)\s+potential\s+at\s+core""", reX)
    for i, line in enumerate(lines[::-1]):
      if regex.search(line) is not None: break
    assert -i + 2 < len(lines), RuntimeError("Could not find average atomic potential in file.")
    regex = compile(r"""(?:\s|\d){8}\s*(-?\d+\.\d+)""")
    result = []
    for line in lines[-i+2:]:
      data = line.split()
      if len(data) == 0: break
      result.extend([m.group(1) for m in regex.finditer(line)])
        
    return array(result, dtype="float64") * eV

  @property 
  @make_cached
  def electronic_dielectric_constant(self):
    """ Electronic contribution to the dielectric constant. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from re import M as multline
    from numpy import array
    regex = r"\s*MACROSCOPIC\s+STATIC\s+DIELECTRIC\s+TENSOR\s*\(including local field effects in DFT\)\s*\n"\
            r"\s*-+\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*-+\s*\n"
    result = self._find_last_OUTCAR(regex, multline)
    assert result is not None, RuntimeError('Could not find dielectric tensor in output.')
    return array([result.group(i) for i in range(1,10)], dtype='float64').reshape((3,3))

  @property 
  @make_cached
  def ionic_dielectric_constant(self):
    """ Ionic contribution to the dielectric constant. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from re import M as multline
    from numpy import array
    regex = r"\s*MACROSCOPIC\s+STATIC\s+DIELECTRIC\s+TENSOR\s+IONIC\s+CONTRIBUTION\s*\n"\
            r"\s*-+\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*-+\s*\n"
    result = self._find_last_OUTCAR(regex, multline)
    assert result is not None, RuntimeError('Could not find dielectric tensor in output.')
    return array([result.group(i) for i in range(1,10)], dtype='float64').reshape((3,3))

  @property 
  def dielectric_constant(self):
    """ Dielectric constant of the material. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    return  self.electronic_dielectric_constant + self.ionic_dielectric_constant

  @property
  @make_cached
  def stresses(self):
    """ Returns total stress at each relaxation step. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import zeros, abs, dot, all, array
    from numpy.linalg import det
    from quantities import eV, J, kbar
    from re import finditer, M 
    if self.isif < 1: return None
    pattern = """\s*Total\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"""\
              """\s*in kB\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"""
    result = []
    with self.__outcar__() as file:
      for regex in finditer(pattern, file.read(), M): 
        stress = zeros((3,3), dtype="float64"), zeros((3,3), dtype="float64") 
        for i in xrange(2): 
          for j in xrange(3): stress[i][j, j] += float(regex.group(i*6+1+j))
          stress[i][0,1] += float(regex.group(i*6+4))
          stress[i][1,0] += float(regex.group(i*6+4))
          stress[i][1,2] += float(regex.group(i*6+4))
          stress[i][2,1] += float(regex.group(i*6+4))
          stress[i][0,2] += float(regex.group(i*6+4))
          stress[i][2,0] += float(regex.group(i*6+4))
        if sum(abs(stress[0].flatten())) > sum(abs(stress[1].flatten())):
          result.append( stress[0] * float(eV.rescale(J) * 1e22\
                         / abs(det(self.structure.cell*self.structure.scale))) * kbar)
        else: result.append(stress[1])
    return array(result)*kbar
  @property
  def stress(self):
    """ Returns final total stress. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    return self.stresses[-1]

  @property
  @make_cached
  def forces(self):
    """ Forces on each atom. """
    if not self.is_dft: raise AttributeError('not a DFT calculation.')
    from numpy import array
    from quantities import angstrom, eV
    from re import compile, finditer, M
    result = []
    pattern = """ *POSITION\s*TOTAL-FORCE\s*\(eV\/Angst\)\s*\n"""\
              """\s*-+\s*\n"""\
              """(?:\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*\n)+"""\
              """\s*-+\s*\n"""\
              """\s*total drift"""
    with self.__outcar__() as file:
      for regex in finditer(pattern, file.read(), M): pass
    return array([u.split()[3:] for u in regex.group(0).split('\n')[2:-2]], dtype="float64")\
           * eV / angstrom



  def __dir__(self):
    result = [ 'ialgo', 'algo', 'is_dft', 'is_gw', 'encut', 'functional',
               'success', 'datetime', 'initial_structure', 'structure',
               'LDAUType', 'HubbardU_NLEP', 'pseudopotential', 'volume',
               'reciprocal_volume', 'density', 'stoichiometry', 'species',
               'isif', 'nsw', 'ismear', 'sigma', 'relaxation', 'smearing',
               'ibrion', 'potim', 'lorbit', 'isym', 'nupdown', 'lmaxmix',
               'istart', 'icharg', 'precision', 'ediff', 'ediffg', 'kpoints',
               'multiplicity', 'ispin', 'name', 'system', 'ionic_charges',
               'valence', 'nelect', 'extraelectron', 'nonscf', 'lwave', 
               'lvtot', 'nelmdl', 'nelmin', 'nbands' ]
    if self.is_dft: 
      result += [ 'energy_sigma0', 'energies_sigma0', 'all_total_energies', 'fermi0K',
                  'halfmetallic', 'cbm', 'vbm', 'total_energies', 'total_energy', 
                  'fermi_energy', 'moment', 'pressures', 'pressure', 'forces', 'stresses',
                  'stress', 'alphabet', 'xc_g0', 'pulay_pressure', 'fft', 'recommended_fft',
                  'partial_charges', 'magnetization', 'eigenvalues', 'occupations', 'electropot', 
                  'electronic_dielectric_constant', 'ionic_dielectric_constant', 
                  'dielectric_constant' ]
    elif self.is_gw:
      result += ['eigenvalues', 'qp_eigenvalues', 'self_energies', 'occupations']
    return result


                 









