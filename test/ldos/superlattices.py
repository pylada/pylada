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

def scale(structure):
  """ Approximate scale of a given structure. """
  n = len(structure.atoms)
  nSi = len([a.type for a in structure.atoms if a.type == 'Si'])
  return float(nSi) / float(n) * structure.scale + float(n - nSi) / float(n) * 5.69

def direction001(n):
  """" Cell and direction for 001. """
  a, b = (n/2, 0) if n % 2 == 0 else (n/2+0.5, 0.5)
  return [ [a, 0, 0], [b, 0.5, 0.5], [0, -0.5, 0.5] ], [1, 0, 0]
def direction011(n):
  """" Cell and direction for 011. """
  return [[1, 0.0, 0.0],[0,0.5,n*0.5], [0,-0.5,n*0.5]], [0, 1, 1]
def direction111(n):
  """" Cell and direction for 111. """
  return [[0.5*n, -0.5, 0.5],[0.5*n,0, -0.5], [0.5,0.5,0]], [1, 1, 1]

def create_sl(path, direction, nmin, nmax, nstep, x=0.5, density=10e0, input='input.py'):
  """ Creates dictionary of superlattices.
  
      :Parameters: 
        path : str 
          Path to output dictionary.
        direction : callable(int)->(3x3, 3x1)
          A callable taking the number of layers on input and returning a tuple
          consisting of a supercell with the correct number of layers and the
          direction of growth of the superlattice. The supercell must be
          have the correct periodicity: two vectors should be parallel to the
          substrate such that layers can actually be defined. Otherwise, there
          is not one-to-one correspondance between atoms and layers: the same
          atom could belong to different layers. A number of examples are given
          in this module: `direction001`, `direction011`, `direction111`.
        nmin : int
          Minimum number of layers.
        nmax : int 
          Maximum number of layers (excluded).
        nstep : int
          Step between layers: [``nmin``, ``nmin``+``nstep``,
          ``nmin``+2*``nstep``, ..., ``nmax``[
        x : float
          Concentration in Si of the superlattice. Should be between 0 and 1.
        density : float
          Kpoint density for escan calculations,
        input : str
          Path to input file containing escan functional.
  """
  from IPython.ipapi import get as get_ipy
  from numpy.linalg import norm, inv
  from pylada.jobs import JobFolder
  from pylada.escan import read_input, exec_input, ReducedKDensity
  from pylada.crystal.binary import zinc_blende
  from pylada.crystal import layer_iterator

  ip = get_ipy()

  input = read_input(input)
  kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

  lattice = zinc_blende()
  lattice.sites[0].type = 'Si', 'Ge'
  lattice.sites[1].type = 'Si', 'Ge'
  lattice.scale = 5.45
  lattice.find_space_group()

  density = 10e0 * max([1e0/norm(u) for u in inv(lattice.cell)])

  rootjobs = ip.user_ns.get('current_jobfolder', JobFolder())
  for n0 in range(nmin, nmax, nstep):
    # create structure
    cell, dir = direction(n0)
    structure = lattice.to_structure(cell) # reduction not necessary if cell function done right.
    N0 = int(len([0 for layer in layer_iterator(structure, dir)]) * x+1e-8)
    for i, layer in enumerate(layer_iterator(structure, dir)):
      for atom in layer: atom.type = 'Si' if i < N0 else 'Ge'
    xp = float(len([0 for atom in structure.atoms if atom.type == 'Si']))
    xp /= float(len(structure.atoms))
    assert abs(xp - x) < 1e-12
    # name and scale.
    structure.name = "{0[0]}{0[1]}{0[2]}/x_{1:0<4.3}/n_{2}".format(dir, x, n0)
    structure.scale = scale(structure)
    # creates job-folder.
    jobfolder = rootjobs / structure.name
    jobfolder.jobparams['structure'] = structure.copy()
    jobfolder.functional = kescan.copy()
    jobfolder.functional.kpoints = ReducedKDensity(density, (0.5, 0.5, 0.5))
    jobfolder.functional.reference = None
    jobfolder.functional.fft_mesh = fftmesh(structure.cell)
    jobfolder.functional.nbstates = int(len(structure.atoms) * 4 * 1.5+0.5)

  if 'current_jobfolder' not in ip.user_ns: ip.user_ns["current_jobfolder"] = rootjobs
  ip.magic("savejobs " + path)
  return
