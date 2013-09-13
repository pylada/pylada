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

""" Tests launching jobs through magic function. """
from pylada.opt import AbstractExtractBase
from pylada.opt.decorators import broadcast_result

def functional(outdir=None, comm=None, **kwargs):
  """ Dummy functional. """
  from os import getcwd
  from os.path import join
  from pylada.opt import Changedir
  if outdir is None: outdir = getcwd()
  comm.barrier()
  with Changedir(outdir, comm) as cwd:
    if comm.is_root:
      with open('out', 'w') as file:
        file.write("comm size : {0}\n".format(comm.size))
        file.write("kargs : {0}\n".format(kwargs))

class Extract(AbstractExtractBase):
  """ Defines success for dummy functional. """
  @property
  @broadcast_result(attr=True, which=0)
  def success(self): 
    """ True if relevant file exists. """
    from os.path import exists, join
    return exists(join(self.directory, 'out'))

functional.Extract = Extract

def create_jobs(n=3):
  """ Returns dictionary with fake jobs. """
  from launch import functional
  from pylada.jobs import JobFolder

  root = JobFolder()
  for i in xrange(n):
    job = root / "results" / str(i)
    job.functional = functional
    job.jobparams["param"] = 'param'

  return root
