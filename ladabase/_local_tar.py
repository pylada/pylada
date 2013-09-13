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

""" Actually performs local transfer 

    As a stand-alone script so we can call from Popen with nohup.
"""
from sys import stderr, exit
from os import fork, getcwd, chmod, remove
from stat import S_IRGRP, S_IWGRP, S_IRUSR, S_IWUSR
from os.path import relpath, dirname
from shutil import move
from pylada.ladabase.ipython import _walk_files, _get_local_push_parser
import tarfile 

# first fork to orphan the process created in the second.
parser = _get_local_push_parser()
parser.add_argument('--tarfile', type=str)
parser.add_argument('--commentfile', type=str)

try: args = parser.parse_args()
except SystemExit: exit(0)

# allows us to exit from Popen call.
try: 
  pid = fork()
  if pid > 0: exit(0)
except OSError as error: 
  print >>stderr, "Fork failed when daemonizing local push: {0.errno} ({0.strerror}.".format(error)
  exit(1)

# now gets list of files.
tarme = tarfile.open(args.tarfile + "_partial" + ".tgz", 'w:gz')
try:
  for rootdir, file in _walk_files(args):
    print file
    tarme.add(file, arcname=relpath(file, dirname(rootdir)))
  tarme.close()
except: 
  remove(args.tarfile + "_partial" + ".tgz")
  remove(args.tarfile + "_partial" + ".comment")
  print >>stderr, "Error while searching and taring output files."
  exit(1)

move(args.tarfile + "_partial" + ".tgz", args.tarfile + ".tgz")
chmod(args.tarfile + ".tgz", S_IRGRP|S_IWGRP|S_IRUSR|S_IWUSR)
chmod(args.tarfile + ".log", S_IRGRP|S_IWGRP|S_IRUSR|S_IWUSR)
chmod(args.tarfile + ".comment", S_IRGRP|S_IWGRP|S_IRUSR|S_IWUSR)
