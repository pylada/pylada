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

def fftmesh(cell, cutoff=20e0):
  """ FFT real-space mesh used by escan. """
  from numpy.linalg import norm
  result = int(norm(cell[:,0]) * cutoff + 5e-1), \
           int(norm(cell[:,1]) * cutoff + 5e-1), \
           int(norm(cell[:,2]) * cutoff + 5e-1)
  # This line makes sure mesh is even. Not sure why...
  result = result[0] + result[0] % 2, \
           result[1] + result[1] % 2, \
           result[2] + result[2] % 2
  return result

def create_start(path, nall = 3, nrand = 5, nmax=100, density=10e0, input='input.py'):
  """ Creates dictionary with input structures for Si/Ge. 
  
      :Parameters:
        path : str
          Path to output dictionary.
        nall : int
          All structure with ``nall`` (excluded) unit-cells are included in the final dictionary.
        nrand : int
          Structures between ``nall`` (included) and ``nrand`` (excluded) unit-cells are
          also considered for inclusion in the final dictionary. However, only
          ``nmax`` are randomly chosen in the end.
        nmax : int
          Structures between ``nall`` (included) and ``nrand`` (excluded) unit-cells are
          also considered for inclusion in the final dictionary. However, only
          ``nmax`` are randomly chosen in the end.
        density : float
          Kpoint density for escan calculations,
        input : str
          Path to input file containing escan functional.
       
      Creates a job-dictionary with a number of structures sampled from an
      exhaustive list of structures to evaluate using escan.
  """
  from random import shuffle
  from itertools import chain
  from IPython.ipapi import get as get_ipy
  from numpy.linalg import norm, inv
  from pylada.enumeration import Enum
  from pylada.crystal.binary import zinc_blende
  from pylada.jobs import JobFolder
  from pylada.escan import read_input, exec_input, ReducedKDensity
  from pylada.crystal.gruber import Reduction

  input = read_input(input)
  kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

  enum = Enum(zinc_blende())
  enum.sites[0].type = 'Si', 'Ge'
  enum.sites[1].type = 'Si', 'Ge'
  enum.scale = 5.45
  enum.find_space_group()
  density = density * max([1e0/norm(u) for u in inv(enum.cell * enum.scale).T])

  strs = [u for  n in range(nall, nrand) for u in enum.xn(n)]
  shuffle(strs)
  strs = [enum.as_structure(*u) for u in strs[:nmax]]
  alls = [structure for n in range(nall) for structure in enum.structures(n)]

  jobs = JobFolder()
  for i, structure in enumerate(chain(alls, strs)):
    structure.name = str(i)
    nSi = len([a.type for a in structure.atoms if a.type == 'Si'])
    structure.scale = float(nSi) / float(n) * enum.scale + float(n - nSi) / float(n) * 5.69

    jobfolder = jobs / structure.name
    jobfolder.jobparams['structure'] = structure.copy()
    jobfolder.structure.cell = Reduction()(jobfolder.structure.cell)
    jobfolder.functional = kescan.copy()
    jobfolder.functional.kpoints = ReducedKDensity(density, (0.5, 0.5, 0.5))
    jobfolder.functional.reference = None
    jobfolder.functional.fft_mesh = fftmesh(structure.cell)
    jobfolder.functional.nbstates = int(len(structure.atoms) * 4 * 1.5+0.5)

  ip = get_ipy()
  ip.user_ns["current_jobfolder"] = jobfolder.root
  ip.magic("savejobs " + path)
  return
