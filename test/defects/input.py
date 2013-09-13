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

""" Input script for the calculation of Point Defects. """
from pylada.crystal import A2BX4

lattice = A2BX4.b5()
""" Back-bone lattice. """
# changes species in lattice.
for site in lattice.sites:
  site.type = {"A": "Rh", "B": "Zn", "X": "O"}[site.type[0]]
lattice.scale = 8.506

supercell = array([[-1,  1,  1],\
                   [ 1, -1,  1],\
                   [ 1,  1, -1]], dtype="float64" )
""" Supercell of defect structures. """

vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
vasp.lorbit     = 10
vasp.npar       = 2
vasp.lplane     = True
vasp.addgrid    = True
vasp.relaxation = "ionic"
vasp.set_smearing   = "metal", 0.01
vasp.set_symmetries = 1e-3

#                Symbol, directory of POTCAR, U parameters, max/min oxidation state, is magnetic
vasp.add_specie = "Fe", "pseudos/Fe", U("liechtenstein", "d", 2.5), 3
vasp.add_specie = "Al", "pseudos/Al", None,  2
vasp.add_specie =  "O",  "pseudos/O", None, -2
vasp.add_specie = "Li", "pseudos/Li", None,  1
vasp.species["Fe"].moment = [5e0, 1e0]


first_trial = { "kpoints": "\n0\nAuto\n1", "encut": 0.9 }
""" parameter to override during first relaxation step. """
relaxation_dof = "volume ionic cellshape"
""" Degrees of freedom to relax. """
relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)
""" Cell shape relaxation algorithm. """


def scale(structure):
  """ Returns *guessed* scale (eg volume^(0.33)) for a given structure. """
  from numpy.linalg import det
  if "O" in [atom.type for atom in structure.atoms]:    spvol = 8.5**3/4e0
  elif "Se" in [atom.type for atom in structure.atoms]: spvol = 9.5**3/4e0
  elif "Te" in [atom.type for atom in structure.atoms]: spvol = 10.5**3/4e0
  else: raise ValueError("Neither O, nor Se, nor Te atom found.")

  nfu = float(len(structure.atoms)/7)*0.5 # 0.5 because 2 f.u. in spinel unit-cell.
  vol = det(structure.cell)
  return (nfu * spvol / vol)**(1e0/3e0) 

materials = ["Al2MgO4"]
""" Materials to compute. """
nbantiferro = 3
""" Number of random anti-ferro trials. """

lattices = [ A2BX4.b4() ]
""" Lattice for which to create structures. """

# The following sets-up which point defects to do.
# point_defects = {"Mn": ["Ga", None], "Ga": ["Mn", None], "O": [None]}
point_defects = {'cations1': 'cations'}
# The following sets-up which which interstitials to examine.
# Li interstitial on 16c (0,0,0)
# Li interstitial on 32e (x,x,x) with x = 0.75.
# These are *FRACTIONAL* coordinates of the LATTICE unit-cell.
# point_defects[None] =\
#      [ "Li", (0,0,0), "16c" ], \
#      [ "Li", (0.75,0.75,0.75), "32e_0.75" ]
# This is a sequence. If there is only one interstitial, use the following line:
# point_defects[None] = [  [ "Li", (0,0,0), "16c" ]  ]

