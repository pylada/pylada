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

""" Submodule defining the *third wave* to compute point-defects. """
__docformat__ = "restructuredtext en"
__all__ = ['pointdefect_wave']

coord_tolerance = 0.25
def pointdefect_wave(path=None, inputpath=None, **kwargs):
  """ Creates point-defect wave using ground-state job-dictionary. 

      :Parameters:
        path : str or None
          Path where the modified job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. If None, will use the
          current job-dictionary path.
        inputpath : str or None
          Path to an input file. If not present, then no input file is read and
          all parameters are taken from the non-magnetic wave.
        kwargs
          Any keyword/value pair to take precedence over anything in the input file.

      Creates a point-defect wave from the materials computed in the
      magnetic and non-magnetic waves. Usage is fairly simple. If the pickle
      for the magnetic/non-magnetic wave is called ``magnetic_wave``, then one
      need only open it and call the ``pointdefect_wave``.

      >>> explore all magnetic_wave 
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered

      The above will add point-defect calculations for all meterials and
      lattices of the ``magnetic_wave`` job-dictionary and save it (to the same 
      path unless an argument is provided to ``magnetic_wave``). Note that
      changing the location of the current job-dictionary has no effect. It
      would be possible but sounds too error prone:

      >>> explore all magnetic_wave 
      >>> goto /Fe2AlO4 # has no effect on what point-defects are added.
      >>> goto /Al2FeO4 # has no effect on what point-defects are added.
      >>> import test
      >>> # creates point-defects for both Fe2AlO4 and Al2FeO4: location does not matter.
      >>> test.magnetic_wave() 
      >>> launch scattered

      Point-defects calculations are added to
      a material if and only if all existing magnetic/non-magnetic jobs for
      that material have completed successfully. Furthermore, only untagged
      materials are accepted. Hence, to disable Fe2AlO4 lattices from having
      point-defects added to it, one can simply tag it:

      >>> explore all magnetic_wave 
      >>> jobparams["/Fe2AlO4"] = 'off'
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered

      Similarly to make sure point-defect calculations are *not* created for
      the b5 lattice of the Fe2AlO4 material, one could tag it as follows:

      >>> explore all magnetic_wave 
      >>> current_jobfolder["/Fe2AlO4/b5"].tag()
      >>> import test
      >>> test.magnetic_wave()
      >>> launch scattered
  """
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  from IPython.ipapi import get as get_ipy
  from numpy import array, sum, abs
  from pylada.jobs import JobFolder
  from pylada.vasp import read_input
  from pylada.opt import Input
  from pylada.crystal import defects as ptd

  # Loads job-folder and path as requested. 
  ip = get_ipy()
  if "current_jobfolder" not in ip.user_ns: 
    print "No current job-dictionary." 
    return
  jobfolder = ip.user_ns["current_jobfolder"].root
  if path is None:
    if "current_jobfolder_path" not in ip.user_ns:
      print "No known path for current dictionary and no path specified on input."
      return
    path = ip.user_ns["current_jobfolder_path"]
  basedir = dirname(path)

  # create input dictionary. First reads non-magnetic input, then magnetic
  # input, then kwargs. Order of precedence is exact opposite.
  input = Input()
  if hasattr(jobfolder, "nonmaginput"):
    with NamedTemporaryFile() as file: 
      file.write(jobfolder.nonmaginput)
      file.flush()
      input.update(read_input(file.name))
  if hasattr(jobfolder, "maginput"):
    with NamedTemporaryFile() as file: 
      file.write(jobfolder.nonmaginput)
      file.flush()
      input.update(read_input(file.name))
  if inputpath is not None:
    input.update(read_input(inputpath))
    with open(inputpath, "r") as file: jobfolder.maginput = file.read()
  input.update(kwargs)

  # saves inputfile to jobfolderioanry if needed.
  if inputpath is not None:
    input.update(read_input(inputpath))
    with open(inputpath, "r") as file: jobfolder.pointdefectinput = file.read()
  # saves current script tof file.
  with open(__file__, "r") as file: jobfolder.pointdefectscript = file.read()

  if hasattr(input, 'coord_tolerance'): coord_tolerance = input.coord_tolerance
  
  nb_new_jobs = 0
  # loops over completed structural jobs.
  for name in magnetic_groundstates():
    # gets the ground-states job-dictionary.
    groundstate = jobfolder[name]
    # checks that the lattice and material are not tagged. 
    if groundstate[".."].is_tagged: continue
    if groundstate["../.."].is_tagged: continue
    # extracts the structure from it.
    superstructure, lattice = create_superstructure(groundstate, input)
    # extracts material.
    material = groundstate.material
    print 'Working on ', groundstate.name, groundstate.name.split('/')[1]
    # extracts description of species.
    species = groundstate.functional.vasp.species

    # loop over substitutees.
    for structure, defect, B in ptd.iterdefects(superstructure, lattice, input.point_defects, tolerance=coord_tolerance):
      # loop over oxidations states.
      for nb_extrae, oxname in charged_states(input.vasp.species, material, defect.type, B):
        
        # creates list of moments. 
        new_moments = deduce_moment(defect.type, species) 
        if len(new_moments) > 1: 
          moments = [ (min(new_moments), "/ls"), (max(new_moments), "/hs") ]
        else:
          moments = [ (max(new_moments), "") ]
        # loop  over moments.
        for moment, suffix in moments:
          name =  "PointDefects/{1}/{2}{0}".format(suffix, structure.name, oxname)
          
          # checks if job already exists. Does not change job if it exists!
          if name in groundstate[".."]: continue

          # creates new job.
          jobfolder = groundstate["../"] / name
          jobfolder.functional = input.relaxer
          jobfolder.jobparams  = groundstate.jobparams.copy()
          jobfolder.jobparams["structure"] = structure.deepcopy()
          jobfolder.jobparams["nelect"] = nb_extrae
          jobfolder.jobparams["relaxation"] = "ionic"
          jobfolder.jobparams["ispin"] = 2
          jobfolder.jobparams["set_symmetries"] = "off"
          jobfolder.lattice  = lattice
          jobfolder.material = material
          jobfolder.defect   = defect
          # adds, modifies, or remove moment depending on defect type.
          if hasattr(superstructure, "magmom") or abs(moment) > 1e-12: 
            jstruct = jobfolder.jobparams["structure"]
            # construct initial magmom
            if hasattr(superstructure, "magmom"):
              jstruct.magmom = [u for u in superstructure.magmom]
            else: 
              jstruct.magmom = [0 for u in superstructure.atoms]
            # now modifies according to structure.
            if B is None or B.lower() == 'none': # interstitial:
              jstruct.magmom.append(moment)
            elif defect.type is None or defect.type.lower() == 'none':
              vacancy_moment(input.vasp.species, jstruct, defect, B, nb_extrae)
            else: 
              jstruct.magmom[defect.index] = moment
            # only keep moment if there are moments. 
            if sum(abs(jstruct.magmom)) < 1e-12 * float(len(jstruct.atoms)): del jstruct.magmom

          nb_new_jobs += 1

  # now saves new job folder
  print "Created {0} new jobs.".format(nb_new_jobs)
  if nb_new_jobs == 0: return
  ip.user_ns["current_jobfolder"] = jobfolder.root
  ip.magic("savejobs " + path)
            
          

def create_superstructure(groundstate, input):
  """ Creates a superstructure from existing structure. """
  from os.path import dirname, join
  from operator import itemgetter
  from numpy import dot
  from IPython.ipapi import get as get_ipy
  from pylada.crystal import fill_structure, vasp_ordered

  # sanity checks,
  assert "structure" in groundstate.jobparams,\
         ValueError("Could not find structure in ground-state job-dictionary.")
  assert hasattr(groundstate.functional, "Extract"),\
         ValueError("Could not find extraction class in ground-state job-dictionary.")
  
  ip = get_ipy()
  assert "current_jobfolder_path" in ip.user_ns,\
         RuntimeError("Could not find path for current job-dictionary.")
  rootdir = dirname(ip.user_ns["current_jobfolder_path"])
  # gets original lattice from job-dictionary.
  orig_structure = groundstate.jobparams["structure"]
  
  # Extracts computed lattice from ground state calculation.
  extract = groundstate.functional.Extract( join(rootdir, groundstate.name[1:]) )
  assert extract.success, RuntimeError("Ground-state computation was not successful.")
  lattice = extract.structure.to_lattice()

  # creates superstructure.
  cell = dot(lattice.cell, input.supercell)
  result = fill_structure(cell, lattice)
  assert len(result.atoms) != len(lattice.sites), \
         ValueError("Superstructure as large as lattice. Disable this line if that's ok.")

  # adds magnetic moment if necessary.
  if hasattr(orig_structure, "magmom"):
    assert extract.magnetization.shape[0] == len(lattice.sites),\
           RuntimeError("Could not find magnetization in ground-state's OUTCAR.")
    mlat = lattice.deepcopy()
    for atom, m in zip(mlat.sites, extract.magnetization[:,-1]):
      if abs(m) < 0.1: atom.type = '0'
      elif m < 0e0: atom.type = str(int(m-1))
      else: atom.type = str(int(m+1))
    moments = fill_structure(cell, mlat)
    result.magmom = [ int(i.type) for i in moments.atoms ]

  return vasp_ordered(result, attributes=["magmom"]), lattice

def magnetic_groundstates():
  """ Yields name of magnetic-groundstates from current job-dictionary.

      A set of magnetic-states for the same lattice and materials is defined by
      all jobs residing in the parent directory of .*/.*/non-magnetic, other
      than PointDefects. 

      All jobs within a set of magnetic-states must be
      finished. Otherwise, that particular combination of material + lattice is
      not considered. 

      This yields the fully qualified job-name of each lowest energy magnetic
      ground-state within the current job-dictionary.
  """
  from operator import itemgetter
  from pylada.jobs import MassExtract
  collect = MassExtract(naked_end=False)
  # loops over untagged non-magnetic structural jobs.
  for nonmag in collect.children:
    # check for success of all jobs (except for Point-defects).
    success = [u[1] for u in nonmag.success.items() if u[0].find("PointDefects") == -1]
    if not all(success): continue
    # checks for lowest energy structure.
    energies = [u for u in nonmag.total_energy.items() if u[0].find("PointDefects") == -1]
    energies = sorted(energies, key=itemgetter(1))
    yield energies[0][0]

def deduce_moment(type, species):
  """ Returns moment.

      This is a helper function which all atomic species the same with respect
      to the attribute ``moment``. If specie has no ``moment`` attribute,
      returns ``[0]``. If it exists and is a scalar, returns ``[moment]``. And
      if already is a list, returns as is.
  """
  if type is None or type.lower() == 'none': return [0]
  if not isinstance(type, str): type = type[0]
  if not hasattr(species[type], "moment"): return [0]
  if not hasattr(species[type].moment, "__iter__"):
    return [species[type].moment]
  return species[type].moment

def oxnumber(species, structure, pos, type):
  """ Returns oxidation number, depending on environment. """
  from pylada.crystal.defects import coordination_number 
  if type not in species: return 0
  if not hasattr(species[type], 'oxidation'): return 0
  result = species[type].oxidation
  if not hasattr(result, '__iter__'): return result
  c = coordination_number(structure, pos, tolerance=coord_tolerance)
  assert c in [4, 6], RuntimeError('Unexpected coordination environment for {0}: {1}'.format(type, c))
  if c == 4: return result[0]
  if c == 6: return result[1]

def charged_states(species, material, Atype, Btype):
  """ List of charged states. """
  from re import compile
  from pylada.crystal.defects import charged_states
  mat_regex = compile('([A-Z][a-z]?)2([A-Z][a-z]?)O4')
  found = mat_regex.match(material)
  A, B = found.group(1), found.group(2)
  Aox = species[Atype].oxidation if Atype is not None and Atype.lower() != 'none' else 0
  Box = species[Btype].oxidation if Btype is not None and Btype.lower() != 'none' else 0
  if hasattr(Aox, '__iter__'): Aox = Aox[1 if Atype == A else 0]
  if hasattr(Box, '__iter__'): Box = Box[1 if Btype == A else 0]

  diff = Box - Aox
  if diff < 0: diff -= 1 
  else: diff += 1 
  sign = 1 if diff > 0 else -1
  for c in range(0, diff, sign):
    if c == 0: yield 0, 'charge_neutral'
    else:      yield c, 'charge_{0}'.format(-c)

def vacancy_moment(species, structure, defect, atomtype, extrae):
  """ Creates/modifies moments surrounding a vacancy. """
  from pylada.crystal.defects import first_shell
  if len(structure.magmom) != len(structure.atoms): structure.magmom.pop(defect.index)
  max_ox = oxnumber(species, structure, defect.pos, atomtype)
  extra_moment = int(abs(max_ox - extrae)+0.001)
  if extra_moment == 0: return
  coordination_shell = first_shell(structure, defect.pos, tolerance=coord_tolerance)
  for n in coordination_shell:
    sign = 1 if structure.magmom[n.index] >= 0e0 else -1
    structure.magmom[n.index] += sign * extra_moment / float(len(coordination_shell))
