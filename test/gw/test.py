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

from pylada.crystal.binary import zinc_blende
from pylada.mpi import world
from pylada.vasp import GWVasp

lattice = zinc_blende()
lattice.scale = 5.45
for site in lattice.sites: site.type = 'Si'

functional = GWVasp()
""" VASP functional """
functional.kpoints      = "Automatic generation\n0\nMonkhorst\n2 2 2\n0 0 0"
functional.precision    = "accurate"
functional.ediff        = 1e-5
functional.encut        = 0.8
functional.lorbit       = 10
functional.npar         = 2
functional.lplane       = True
functional.addgrid      = True
functional.set_smearing = "metal", 0.01
functional.relaxation   = "static"
functional.nbands       = 20
functional.vasp_library = "libvasp-5.2.11.so"

functional.add_specie = "Si", "pseudos/Si"

result = functional(lattice.to_structure(), outdir="results", comm=world)
