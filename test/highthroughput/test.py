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

""" High-Thoughput of A2BO4 structures. """
__docformat__ = "restructuredtext en"
__all__ = ['nonmagnetic_wave', 'magnetic_wave', 'calcGW']

def nonmagnetic_wave(path, inputpath="input.py", **kwargs):
  """ Jobs to explore possible ground-states. 
  
      :Parameters:
        path 
          Path where the job-folder will be saved. Calculations will be
          performed in the parent directory of this file. Calculations will be
          performed in same directory as this file.
        inputpath
          Path to an input file. Defaults to input.py. 
        kwargs
          Any keyword/value pair to take precedence over anything in the input
          file.

      Creates a high-throughput job-folder to compute the non-magnetic
      ground-state of a host-material.  The new job-folder is loaded into
      memory automatically. No need to call explore. It is save to the path
      provided on input.
  """
  import re
  from IPython.core.interactiveshell import InteractiveShell
  from copy import deepcopy
  from pylada.vasp import read_input
  from pylada.jobfolder import JobFolder
  from pylada import interactive

  from pylada.misc import setBugLev
  setBugLev(0)   # set global debug level
  from pylada.misc import bugLev

  # reads input.
  input = read_input(inputpath)
  print '  test/hi/test nonmag: inputpath: ', inputpath
  print '  test/hi/test nonmag: input.vasp.relaxation: ', input.vasp.relaxation
  input.update(kwargs)


#  """ Materials to compute. """
#  materials = [ "Al2MgO4"]
#
#  """ Number of random anti-ferro trials. """
#  nbantiferro = 8
#  nbrandom    = 3
#  do_ferro    = False
#  do_antiferro = False
#
#  from pylada.crystal import A2BX4
#  lattices = [A2BX4.b5(), A2BX4.b21()]
#
#  mlen = len( materials)
#  llen = len( lattices)
#  matLatPairs = (mlen * llen) * [None]
#  print "  test/hi/test nonmag: mlen: ", mlen
#  print "  test/hi/test nonmag: llen: ", llen
#  print "  test/hi/test nonmag: pairs len: ", len(matLatPairs)
#
#  kk = 0
#  for mat in materials:
#    print "  test/hi/test nonmag: mat: ", mat
#    for lat in lattices:
#      print "    test/hi/test nonmag: lat: ", lat
#      matLatPairs[kk] = (mat, lat,)
#      kk += 1
#
#  print "test/hi/test nonmag: mats len: %d  lats len: %d  matLatPairs len: %d" \
#    % (mlen, llen, len( matLatPairs),)









  # Job dictionary.
  jobfolder = JobFolder()

  # loop over material-lattice pairs.
  for (material,lattice) in input.matLatPairs:
    if bugLev >= 1:
      print '  test/hi/test nonmag: start material: ', material
      print '  test/hi/test nonmag: start lattice: ', lattice
      print ''
    if bugLev >= 5:
      print '  test/hi/test nonmag: ========== species =========='
      skeys = input.vasp.species.keys()
      skeys.sort()
      for skey in skeys:
        print '    test/hi/test nonmag: species[%s]: %s' \
          % (skey, input.vasp.species[skey], )
      print ''

    if getattr( input, 'useInputCif', None) == True:
      structure = deepcopy(lattice)

    else:
      # Check material
      regex = "([A-Z][a-z]?)2([A-Z][a-z]?)([A-Z][a-z]?)4"
      match = re.match( regex, material)
      if match == None:
        raise RuntimeError("Incorrect material: \"%s\"" % (material,))

      # Checks species are known to vasp functional
      for i in range(1, 4):
        assert match.group(i) in input.vasp.species,\
          RuntimeError("No pseudo-potential defined for {0}.".format(
            match.group(i)))
      # actually creates dictionary.
      species_dict = {"A": match.group(1), "B": match.group(2),
        "X": match.group(3)}
      #print '  test/hi/test nonmag: species_dict: ', species_dict
      
      # Check lattice name
      if len(getattr(lattice, 'name', '').strip()) == 0:
        raise ValueError("Lattice has no name.")


      # creates a structure.
      structure = deepcopy(lattice)
      # changes atomic species, from ABX to real element names.
      # from: Atom(0.5, 0.5, 0.5, 'A')
      # to:   Atom(0.5, 0.5, 0.5, 'Al')
      for atom in structure:  atom.type  = species_dict[atom.type]
      # assigns it a name.
      structure.name = "{0} in {1}, spin-unpolarized.".format(
        material, lattice.name)

      # gets its scale.
      structure.scale = input.scale(structure)




    # job folder for this lattice.
    lat_jobfolder = jobfolder / material 

    job = lat_jobfolder / lattice.name / "non-magnetic"
    job.functional = input.vasp
    job.params["structure"] = structure
    job.params["ispin"] = 1
    # saves some stuff for future reference.
    job.material = material
    job.lattice  = lattice
    #print '    test/hi/test nonmag: job: ', job
    #print '    test/hi/test nonmag: === job.functional ===\n%s\n=== end functional === ' % (job.functional,)


  interactive.jobfolder = jobfolder
  InteractiveShell.instance().magic("savefolders " + path)








def magnetic_wave( path=None, inputpath='input.py', **kwargs):
  """ Creates magnetic wave for current job-folder.

      :Parameters:
        path : str or None
          Path where the modified job-folder will be saved. Calculations will be
          performed in the parent directory of this file. If None, will use the
          current job-folder path.
        inputpath : str or None
          Path to an input file. If not present, then no input file is read and
          all parameters are taken from the non-magnetic wave.
        kwargs
          Any keyword/value pair to take precedence over anything in the input file.

      Creates magnetic wave from pre-existing non-magnetic wave. If no input
      file is given on input, then all parameters are obtained from the
      corresponding non-magnetic wave. 

      The new job-folder is loaded into memory automatically. No need to
      call explore. It is save to the path provided on input (or to the current
      job-folder path ifnot provided).  It will contain magnetic and
      non-magnetic calculations both. Pre-existing magnetic
      calculations will *not* be overwritten. However, additional anti-ferro
      configurations can be calculated by giving a large enough ``nbantiferro``
      on input.
  """
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  from IPython.core.interactiveshell import InteractiveShell
  from pylada.vasp import read_input
  from pylada.jobfolder import JobFolder
  from pylada import interactive
  import pickle, random

  from pylada.misc import setBugLev
  setBugLev(0)   # set global debug level
  from pylada.misc import bugLev

  with open(path) as fin:
    jobfolder = pickle.load( fin)
  basedir = dirname( path)
  jobfolder_path = join( basedir, jobfolder.name)

  # Loads job-folder and path as requested. 
  if jobfolder is None: 
    print "No current job-folder."
    return
  if jobfolder_path is None: 
    print "No path for current job-folder."
    return

  input = read_input(inputpath)
  print '  test/hi/test mag: inputpath: ', inputpath
  print '  test/hi/test mag: input.vasp.relaxation: ', input.vasp.relaxation

  # will loop over all folders, looking for *successfull* *non-magnetic* calculations. 
  # Only magnetic folders which do NOT exist are added at that point.
  nonmagname = "non-magnetic"
  nb_new_folders = 0
  for name, nonmagjob in jobfolder.iteritems():
    if bugLev >= 1:
      print "test/hi/test mag: name: %s  nonmagjob: %s" % (name, nonmagjob,)
    # avoid tagged folders.
    if nonmagjob.is_tagged: continue
    # avoid other folders (eg magnetic folders).
    basename = normpath("/" + name + "/../")
    if bugLev >= 1:
      print "test/hi/test mag: basename: %s" % (basename,)
    if relpath(name, basename[1:]) != nonmagname: continue
    # check for success and avoid failures.
    extractPath = join(basedir, name)
    # xxx use contcar
    extract = nonmagjob.functional.Extract( extractPath)
    if bugLev >= 1:
      print "test/hi/test mag: extractPath: %s" % (extractPath,)
      print "test/hi/test mag: extract.success: %s" % (extract.success,)
      print "test/hi/test mag: extract: %s" % (extract,)
      print "test/hi/test mag: dir(extract): %s" % (dir(extract),)
      print "test/hi/test mag: extract.functional: %s" % (extract.functional,)
    if not extract.success: continue

    if bugLev >= 1:
      print 'test/hi/test mag ========== A'
      print "test/hi/test mag: nonmag structure: %s" % (nonmagjob.structure,)
      print 'test/hi/test mag ========== B'
      print "test/hi/test mag: extract.structure: %s" % (extract.structure,)
      print 'test/hi/test mag ========== C'
      print "test/hi/test mag: extract.species: %s" % (extract.species,)
      print 'test/hi/test mag ========== D'

    if not is_magnetic_system(extract.structure, extract.functional.species):
      continue

    # loads lattice and material from non-magnetic job.
    material = nonmagjob.material
    lattice = nonmagjob.lattice
    if bugLev >= 1:
      print "test/hi/test mag: material: %s" % (material,)
      print "test/hi/test mag: lattice: %s" % (lattice,)

    # figures out whether we have both high and low spins. 
    if bugLev >= 1:
      print "test/hi/test mag: extract structure: %s" % (extract.structure,)
      print "test/hi/test mag: extract species: %s" \
        % (extract.functional.species,)
    if has_high_and_low(extract.structure, extract.functional.species):
          hnl = [(min, "ls-"), (max, "hs-")]
    else: hnl = [(min, "")] 
    # now loops over moments.
    for func, prefix in hnl: 
      if bugLev >= 1:
        print "test/hi/test mag: func: %s  prefix: %s" % (func, prefix,)
        # Either func=min,prefix="ls-", or func=max,prefix="hs-"
      # Now tries and creates high-spin ferro folders
      # if it does not already exist.
      jobname = normpath("{0}/{1}ferro".format(basename, prefix))
      structure, magmom = ferro(
        extract.structure, extract.functional.species, func)

      if bugLev >= 1:
        print "test/hi/test mag: structure: %s" % (structure,)
        print "test/hi/test mag: jobfolder: %s" % (jobfolder,)
        print "test/hi/test mag: jobname: %s" % (jobname,)
        print "test/hi/test mag: magmom: %s" % (magmom,)
        print "test/hi/test mag: input.do_ferro: %s" % (input.do_ferro,)
      if magmom and jobname not in jobfolder and input.do_ferro:
        structure.name = "{0} in {1}, {2}ferro."\
                         .format(material, lattice.name, prefix)
        job = jobfolder / jobname
        if bugLev >= 1:
          print "test/hi/test mag: structure.name: %s" % (structure.name,)
          print "test/hi/test mag: job: %s" % (job,)
        job.functional = input.vasp
        job.params["structure"] = structure.copy()
        job.params["magmom"] = True
        job.params["ispin"] =  2
        # saves some stuff for future reference.
        job.material = material
        job.lattice  = lattice
        print "test/hi/test mag: created ferro jobname: %s" % (jobname,)
        nb_new_folders += 1

      # Now tries and creates anti-ferro-lattices folders
      # if it does not already exist.
      structure, magmom = species_antiferro(
        extract.structure, extract.functional.species, func) 

      jobname = normpath("{0}/{1}anti-ferro-0".format(basename, prefix))
      if magmom and jobname not in jobfolder and input.do_antiferro:
        structure.name = "{0} in {1}, {2}specie-anti-ferro."\
                         .format(material, lattice.name, prefix)

        job = jobfolder / jobname
        job.functional = input.vasp
        job.params["structure"] = structure.copy()
        job.params["magmom"] = True
        job.params["ispin"] =  2
        # saves some stuff for future reference.
        job.material = material
        job.lattice  = lattice
        print "test/hi/test mag: created anti ferro jobname: %s" % (jobname,)
        nb_new_folders += 1

      # Random anti-ferro, or all possible anti-ferro.
      magIxs = []         # indices of magnetic atoms
      for ii in range( len( extract.structure)):
        atom = extract.structure[ii]
        if len(deduce_moment(atom, extract.functional.species)) > 1:
          magIxs.append( ii)
      numMagAtom = len( magIxs)

      # Total num tests = 2**numMagAtom, but each one is counted
      # twice, since ddd == uuu, and ddu == uud, etc.,
      # where d = spin down, u = spin up.
      # So the num unique tests = 2**numMagAtom / 2
      numUnique = 2 ** (numMagAtom - 1)

      # If the user asked for at least half the total num unique tests,
      # we may as well do them all
      print "test/hi/test mag: numMagAtom: %d  numUnique: %d  input.nbantiferro: %d" \
        % ( numMagAtom, numUnique, input.nbantiferro,)
      if input.nbantiferro >= numUnique / 2:
        # Do all possible combinations
        print 'test/hi/test mag: Do all possible antiferro combinations'
        for itest in range( numUnique):
          jobname = runAntiFerroTest(
            bugLev, jobfolder, material, lattice,
            inputpath, input, magIxs,
            extract, basename, func, prefix, itest, itest)
          print "test/hi/test mag: created all anti ferro jobname: %s" \
            % (jobname,)
          nb_new_folders += 1

      else:
        # Only run a few random combinations
        print 'test/hi/test mag: Do a random selection of antiferro combinations'
        doneMap = {}
        ix = 0
        while True:
          itest = random.randint( 0, numUnique - 1)
          if itest not in doneMap:
            doneMap[itest] = True
            jobname = runAntiFerroTest(
              bugLev, jobfolder, material, lattice,
              inputpath, input, magIxs,
              extract, basename, func, prefix, ix, itest)
            print "test/hi/test mag: created random anti ferro jobname: %s" \
              % (jobname,)
            nb_new_folders += 1
            ix += 1
            if ix >= input.nbantiferro: break

  # now saves new job folder
  print "Created {0} new folders.".format(nb_new_folders)
  if nb_new_folders == 0: return
  interactive.jobfolder = jobfolder.root
  InteractiveShell.instance().magic("savefolders " + path)







def is_magnetic_system(structure, species):
  """ True if system is magnetic. """
  from pylada.misc import bugLev

  # xxx structure has no .atoms
  uatoms = sorted(list(set([a.type for a in structure])))
  if bugLev >= 1:
    print "test/hi/test.is_mag: uatoms: %s" % (uatoms,)
    print "test/hi/test.is_mag: species.items(): %s" % (species.items(),)

  # for each u==Specie in species.items where name is in uatoms ...
  for u in [u for name, u in species.items() if name in uatoms]:
    if bugLev >= 1:
      print "test/hi/test.is_mag: u: %s" % (u,)
      print "test/hi/test.is_mag: hasattr( u, \"moment\"): %s" \
        % (hasattr( u, "momemt"),)
    if not hasattr(u, "moment"): continue

    if bugLev >= 1:
      print "test/hi/test.is_mag: hasattr( u.moment, \"_itier__\"): %s" \
        % (hasattr( u.moment, "__iter__"),)
    if not hasattr(u.moment, "__iter__"): 
      if abs(u.moment) > 1e-12: return True
      continue
    for a in u.moment:
      if abs(a) > 1e-12: return True
    
  return False



def has_high_and_low(structure, species):
  """ True if some species have both high and low spins. """
  for atom in structure:
    if len(deduce_moment(atom, species)) > 1: return True
  return False



def deduce_moment(atom, species):
  """ Returns moment.

      This is a helper function which all atomic species the same with respect
      to the attribute ``moment``. If specie has no ``moment`` attribute,
      returns ``[0]``. If it exists and is a scalar, returns ``[moment]``. And
      if already is a list, returns as is.
  """
  if not hasattr(species[atom.type], "moment"): return [0]
  if not hasattr(species[atom.type].moment, "__iter__"):
    return [species[atom.type].moment]
  return species[atom.type].moment



def ferro(structure, species, func=min):
  """ Returns magmom VASP flag for low-spin ferromagnetic order. """
  result = structure.copy()
  for atom in result: atom.magmom = func(deduce_moment(atom, species))
  return result, True



def species_antiferro(structure, species, func=min):
  """ Low spin anti ferro order with each cation specie in a different direction. """
  from numpy import math, abs
  # checks whether lattice-sites are magnetic or not.
  result = structure.copy()
  # signs = map: atomSym -> 1 if not mag, -1 if mag
  signs = {}
  for atom in result:
    m = func(deduce_moment(atom, species))
    signs[atom.type] = 1 if abs(m) < 1e-12 else -1
  if len([u[1] for u in signs.items() if u[1] == -1]) < 2: return None, False
  # makes alternating sign list.
  # All atoms of a given specie get the same sign.
  dummy = 1
  for k in sorted(signs.keys()):
    if signs[k] != -1: continue
    signs[k] = 1 if dummy == 1 else -1
    dummy = -1 * dummy

  # adds magmom tag to structure's atoms.
  s = 0e0
  for atom in result:
    atom.magmom = float(signs[atom.type]) * func(deduce_moment(atom, species))
    s += abs(atom.magmom)
  return result, s > 1e-6 





# Given itest use it's binary representation to set the spins.
# For example if numMagAtom = 4 and itest = 5,
# we set the spins to 0101.
# Returns jobname.

def runAntiFerroTest(
  bugLev, jobfolder, material, lattice,
  inputpath, input, magIxs,
  extract, basename, func, prefix, ix, itest):

  from os.path import normpath

  structure = extract.structure.copy()
  numMagAtom = len( magIxs)
  if bugLev >= 1:
    print 'runAntiFerroTest: material: %s  ix: %d  itest: %d' \
      % (material, ix, itest,)
    print 'runAntiFerroTest: numMagAtom: %d  magIxs: %s' \
      % (numMagAtom, magIxs,)

  jobname = normpath("{0}/{1}anti-ferro-{2}".format(
    basename, prefix, ix))

  # Shift the bits in itest right, examining the low order bit,
  # and make corresponding array: spins[], having elements -1 and 1.
  spins = numMagAtom * [0]
  itmp = itest
  for ii in range( numMagAtom):
    spin = -1.
    if itmp % 2 == 1: spin = 1.
    spins[ii] = spin
    itmp >>= 1

  # Copy the moments[] array into the magnetic atoms.
  for ii in range( numMagAtom):
    spin = spins[ii]
    atom = structure[ magIxs[ii]]
    atom.magmom \
      = spin * func( deduce_moment( atom, extract.functional.species))

  if jobname not in jobfolder:
    structure.name = "{0} in {1}, random anti-ferro.".format(
      material, lattice.name)
    job = jobfolder / jobname
    if inputpath == None: job.functional = nonmagjob.functional
    else: job.functional = input.vasp
    job.params["structure"] = structure.copy()
    job.params["magmom"] = True
    job.params["ispin"] =  2
    # saves some stuff for future reference.
    job.material = material
    job.lattice  = lattice

  return jobname















def calcGW( path=None, inputpath='inputGW.py', **kwargs):
  """ Runs GW calcs on all subfolders.

      :Parameters:
        path : str or None
          Path where the modified job-folder will be saved. Calculations will be
          performed in the parent directory of this file. If None, will use the
          current job-folder path.
        inputpath : str or None
          Path to an input file.
        kwargs
          Any keyword/value pair to take precedence over anything in the input file.
  """
  from tempfile import NamedTemporaryFile
  from os.path import dirname, normpath, relpath, join
  from IPython.core.interactiveshell import InteractiveShell
  from pylada.vasp import read_input
  from pylada.jobfolder import JobFolder
  from pylada import interactive
  import pickle, random

  from pylada.misc import setBugLev
  setBugLev(0)   # set global debug level
  from pylada.misc import bugLev

  with open(path) as fin:
    jobfolder = pickle.load( fin)
  basedir = dirname( path)
  jobfolder_path = join( basedir, jobfolder.name)

  # Loads job-folder and path as requested. 
  if jobfolder is None: 
    print "No current job-folder."
    return
  if jobfolder_path is None: 
    print "No path for current job-folder."
    return

  input = read_input(inputpath)
  print '  test/hi/test gw: inputpath: ', inputpath
  print '  test/hi/test gw: input.vasp.relaxation: ', input.vasp.relaxation

  # will loop over all folders, looking for all successfull calculations. 
  nb_new_folders = 0
  for name, pjob in jobfolder.iteritems():
    if bugLev >= 1:
      print "test/hi/test gw: name: %s  pjob: %s" % (name, pjob,)
    # avoid tagged folders.
    if pjob.is_tagged: continue

    basename = normpath("/" + name + "/../")
    if bugLev >= 1:
      print "test/hi/test gw: basename: %s" % (basename,)

    # check for success and avoid failures.
    extractPath = join(basedir, name)
    extract = pjob.functional.Extract( extractPath)
    if bugLev >= 1:
      print "test/hi/test gw: extractPath: %s" % (extractPath,)
      print "test/hi/test gw: extract.success: %s" % (extract.success,)
      print "test/hi/test gw: extract: %s" % (extract,)
      print "test/hi/test gw: dir(extract): %s" % (dir(extract),)
    if not extract.success: continue

    if bugLev >= 1:
      print 'test/hi/test gw ========== A'
      print "test/hi/test gw: nonmag structure: %s" % (pjob.structure,)
      print 'test/hi/test gw ========== B'
      print "test/hi/test gw: extract.structure: %s" % (extract.structure,)
      print 'test/hi/test gw ========== C'
      print "test/hi/test gw: extract.species: %s" % (extract.species,)
      print "test/hi/test gw: extract functional species: %s" \
        % (extract.functional.species,)
      print 'test/hi/test gw ========== D'
      print "test/hi/test gw: dir(pjob): %s" % (dir(pjob),)
      print 'test/hi/test gw ========== E'

    # loads lattice and material from previous job.
    material = pjob.material
    lattice = pjob.lattice
    if bugLev >= 1:
      print "test/hi/test gw: material: %s" % (material,)
      print "test/hi/test gw: lattice: %s" % (lattice,)

    jobname = pjob.name + '/gwcalc'
    structure = extract.structure.copy()
    structure.name = "{0} in {1}, GW.".format(material, lattice.name)
    if bugLev >= 1:
      print "test/hi/test gw: jobname: %s" % (jobname,)
      print "test/hi/test gw: structure.name: %s" % (structure.name,)

    job = jobfolder / jobname
    if bugLev >= 1:
      print "test/hi/test gw: job: %s" % (job,)

    job.functional = input.vasp

    job.params["structure"] = structure.copy()
    job.params["magmom"] = True
    job.params["ispin"] =  2
    # saves some stuff for future reference.
    job.material = material
    job.lattice  = lattice
    print "test/hi/test gw: created ferro jobname: %s" % (jobname,)
    nb_new_folders += 1

  # now saves new job folder
  print "Created {0} new folders.".format(nb_new_folders)
  if nb_new_folders == 0: return
  interactive.jobfolder = jobfolder.root
  InteractiveShell.instance().magic("savefolders " + path)






