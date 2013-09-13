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

""" Interactive magic functions for CRYSTAL. """
__docformat__ = "restructuredtext en"
__all__ = [] # empty: nothing is actually imported, everython is registered.
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
  from pylada import is_interactive
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

      It can be a bit difficult to determine whether an MPPcrystal job will run
      correctly (e.g. with k-point parallelization), especially if more than
      one job is involved and each job must have the same number of processors.
   
      This magic function determines likely core counts (per job) across a
      job-folder. It can take a single integer argument to filter the processor
      counts to multiples of that integer (e.g. to run on complete
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

@register_line_magic
def crystal_change_spinedit(line):
  """ Assigns spins configuration according to current optimization step. 
  
      Goes to last (successful) step and checks the spin on the structure.
      It then sets edits the
      :py:attr:`~pylada.dftcrystal.functional.scf.atomspin` tag to reproduce the
      computed spin density configuration. 

      The option '--print' will only print the changes without effecting them.
  """
  from operator import attrgetter
  shell = get_ipython()

  # figure out jobparams
  if 'jobparams' not in shell.user_ns:
    print 'No jobfolder loaded in memory.'
    return
  jobparams = shell.user_ns['jobparams']
  collect = shell.user_ns['collect']

  line = line.rstrip().lstrip()
  if '--print' in line: 
    line = line.replace('--print', '').rstrip().lstrip()
    doprint = True
  else: doprint = False
  if '-h' in line.split() or '--help' in line.split():
    print globals()['crystal_change_spinedit'].__doc__
    return
  if len(line) == 0: cutoff = 0e0
  elif line in shell.user_ns: cutoff = float(shell.user_ns[line])
  else: cutoff = float(eval(line, shell.user_ns.copy()))

  for name, params in jobparams.iteritems():
    if params.is_tagged: continue
    if not params.functional.dft.spin: continue
    doup = params.functional.atomspin.up is not None
    dodown = params.functional.atomspin.down is not None
    if not (doup or dodown): continue
    try: details = collect[name].details
    except: 
      print 'bypassing', name
      continue

    if len(details) == 0: continue
    candidates = [n for n, u in details.iteritems() if u.success]
    if len(candidates) == 0: continue
    def cmp(a): return int(a.split('/')[-1])
    extract = details[max(candidates, key=cmp)]


    if doup: 
      atomup = [a for a in extract.structure if a.spin > cutoff]
      atomup = sorted(atomup, key=attrgetter('spin'))[::-1]
      target = len(params.functional.atomspin.up)
      if len(atomup) > target: 
        result = atomup[:target]
        for u in atomup[target:]:
          if abs(u.spin - result[-1].spin) < 1e-4: result.append(u)
        atomup = result
      atomup = sorted(atomup, key=attrgetter('label'))
      if doprint:
        atomspin = [float(u.spin) for u in atomup]
        atomup = [u.label for u in atomup]
        print "{0} (up): {1} -> {2} ({3})"                                     \
              .format(name, params.functional.atomspin.up, atomup, atomspin)
      else: params.functional.atomspin.up = [u.label for u in atomup]
    if dodown: 
      atomdown = [a for a in extract.structure if a.spin < cutoff]
      atomdown = sorted(atomdown, key=attrgetter('spin'))
      target = len(params.functional.atomspin.down)
      if len(atomdown) > target: 
        result = atomdown[:target]
        for u in atomdown[target:]:
          if abs(u.spin - result[-1].spin) < 1e-4: result.append(u)
        atomdown = result
      if doprint:
        atomdown = [u.label for u in atomdown]
        print "{0} (down): {1} -> {2}"                                         \
              .format(name, params.functional.atomspin.down, atomdown)
      else: params.functional.atomspin.down = [u.label for u in atomdown]


@register_line_magic
def crystal_equivlabels(line):
  """ Atomic labels equivalent to those in input. 
  
      Prints the labels which are symmetric equivalent to those given on input.
      
      >>> %crystal_equivlabels 5 6
      5 -> 6, 31, 32
      6 -> 5, 31, 32 

      The above prints the atoms which equivalent to atom 5 in the CRYSTAL_
      structure, and to atom 6.
  """
  from re import split
  from ..crystal.iterator import equivalence as equivalence_iterator
  from ..dftcrystal import Molecule
  from ..interactive import jobfolder as cjf

  if 'structure' not in cjf.params:
    print 'No structure in ', cjf.name
    return

  try: labels = [int(u) for u in split('\s*,?', line) if len(u)]
  except: 
    print "Could not translate input to integers."
    print line
    return

  structure = cjf.structure
  if not isinstance(structure, Molecule):
    print "Structure in current folder is not for CRYSTAL."
    return 
  symops = structure.symmetry_operators
  structure = structure.eval()
  atoms = [ [structure[u].label for u in v]
            for v in equivalence_iterator(structure, symops)]

  for label in labels:
    for group in atoms:
      if label in group: 
        print label, '->', ', '.join(str(u) for u in group if u != label)
        break
  return

