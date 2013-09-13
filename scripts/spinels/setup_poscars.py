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


#! /uhome/mdavezac/usr/bin/python
# 
#PBS -l nodes=2:ppn=2,walltime=48:00:0
#PBS -q Std
#PBS -m n 
#PBS -e err
#PBS -o out

def create_structure_dir(structure, prefix=None):
  from os import makedirs
  from os.path import join, exists

  # create structure name
  nAl = len([ 0 for v in structure.atoms if v.type == "Al" ])
  nMg = len([ 0 for v in structure.atoms if v.type == "Mg" ])
  nO = len([ 0 for v in structure.atoms if v.type == "O" ])
  structure.name = "Al%iMg%iO%i" % ( nAl, nMg, nO )

  # create directory
  dirname = structure.name
  if prefix is not None: dirname = join(prefix, dirname)
  u = 1
  while exists(dirname) :
    dirname = "%s-%i" % (structure.name, u) 
    if prefix is not None: dirname = join(prefix, dirname)
    u += 1

  # create directory, POSCAR, and sendme 
  makedirs(dirname)
  return dirname

def choose_structures( howmany, filename="database", prefix=None ):
  from random import randint
  from pylada import crystal, enumeration
  import database

  N = 0
  for structure in database.read_database(filename, withperms=True): N  += 1

  which = set([])
  while len(which) != howmany: which.add(randint(0,N-1))
  which = sorted(list(which))
  for n, structure in enumerate(database.read_database(filename, withperms=True)):
    if n != which[0]: continue

    which.pop(0)

    dirname = create_structure_dir(structure, prefix)
    structure.name += ", database %s, structure %i" % (filename, n)
    crystal.print_poscar(structure, ("Al", "Mg", "O"), dirname )

    if len(which) == 0: break
    

if __name__ == "__main__":
  choose_structures(35)
