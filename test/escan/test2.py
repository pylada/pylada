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

from os.path import join, exists
from os import makedirs
from shutil import rmtree, copy
from boost.mpi import world
from pylada.escan import call_escan as escan_as_library
from pylada.opt.changedir import Changedir
from pylada.opt import Redirect

# pseudo files
pseudo_files = [ "maskr", "vq.Ge", "vq.Si", "vq.SiGe.Ge", "vq.SiGe.Si", "vwr.pso" ]

# directories where jobs are at.
jobs = ["VBM", "Gamma",  "X",  "L",  "W1", "W2"]
testdir = "test_input"
workdir = "work"
if world.rank == 0 and (not exists(workdir)): makedirs(workdir)

# creates local comm split into N groups.
N = 2
color = world.rank % N
local_comm = world.split(color)

# launch escan for different jobs.
for i, dir in enumerate(jobs):
  # splits job according to color.
  if i % N != color: continue

  # creates working directory with all input files.
  workhere = join(workdir, dir)
  if local_comm.rank == 0:
    if exists(workhere): rmtree(workhere) # deletes everything if already there.
    makedirs(workhere)
    # symlink potential files
    testhere = join(testdir, "pseudos")
    for file in pseudo_files: copy(join(testhere, file), workhere)
    # symlinks input files.
    testhere = join(testdir, dir)
    for file in ["atom_config", "pot_input", "escan_input"]: 
      copy(join(testhere, file), workhere)

  # sync all coms. 
  local_comm.barrier()
  # changes current directory to working directory.
  with Changedir(workhere) as current_dir:
    # Redirects fortran output and error.
    stdout = "stdout"
    stderr = "stderr"
    if local_comm.rank != 0: 
      stdout = "%s.%i" % (stdout, local_comm.rank)
      stderr = "%s.%i" % (stderr, local_comm.rank)
    with Redirect(Redirect.fortran.output, stdout) as stdout:
      with Redirect(Redirect.fortran.error, stderr) as stderr:
        # finally calls escan.
        escan_as_library( local_comm, \
                          atom="atom_config",\
                          pot="pot_input", \
                          escan="escan_input" )
  
