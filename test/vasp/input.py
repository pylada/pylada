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

from quantities import angstrom

structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43)\
                     .add_atom(0,0,0, "Si")\
                     .add_atom(0.25,0.25,0.25, "Si")


vasp = Vasp()
""" VASP functional """
vasp.kpoints    = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
vasp.precision  = "accurate"
vasp.ediff      = 1e-5
vasp.encut      = 1
# vasp.lorbit     = 10
# vasp.npar       = 2
# vasp.lplane     = True
# vasp.addgrid    = True
# vasp.restart_from_contcar = False
vasp.set_smearing   = "metal", 0.01
vasp.relaxation = "volume", 50, 2
# vasp.set_symmetries = "on"

vasp.add_specie = "Si", "pseudos/Si"

# first_trial = { "encut": 0.9 }
# relaxation_dof = "volume ionic cellshape"
# relaxer = RelaxCellShape( vasp, relaxation_dof, first_trial, maxiter=5)
