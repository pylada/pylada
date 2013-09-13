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

# davezac@hopper08:~/SiGe/old/jwluo/SiGe_SL_001/Si6Ge4_on_Ge/vff
from numpy import array
from quantities import angstrom
from pylada.escan import read_input, bandgap
from pylada.crystal import sort_layers
from pylada.crystal.binary import zinc_blende
from pylada.physics import a0

input = read_input("input.py")
cell = [[2.5,0,0],[0.5,0.5,0.5],[0,-0.5,0.5]]
subs = {'A':'Si', 'B':'Ge'}

structure = zinc_blende().to_structure(cell, subs)
structure.scale = float(5.33847592 / 0.5 * a0.rescale(angstrom))
structure = sort_layers(structure, array([2.5,0,0]))
for i, atom in enumerate(structure.atoms):
  atom.type = 'Si' if i % 10 < 6 else 'Ge'

result = bandgap( input.escan, structure, eref=None, 
                  outdir="results/dipoles", references = (-0.2, 0.2),
                  nbstates=4, direction = array([1., 0,0]),
                  fft_mesh = (50,14,14) )

tot = array([0e0] * 3)  / a0 / a0
for e0, e1, u in result.dipole(degeneracy = 5e1 * input.escan.tolerance, attenuate=False):
  tot += (u * u.conjugate()).real
check = array([1.13753549e-03,   5.87932303e-05,   5.87932303e-05]) * 1/a0**2
for a, b in zip(abs(check - tot), abs(check)*1e-2): assert a < max(b, 1e-6), (a, b)
