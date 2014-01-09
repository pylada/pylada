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

import os.path
from glob import iglob
from pylada.crystal.read import icsd_cif_a
from pylada.crystal.read import icsd_cif_b

from pylada.misc import bugLev
from pylada.misc import setTestValidProgram

if bugLev >= 1: print "test/hi/inputCif: entry"

useInputCif = True          # flag used by test.py

vasp = Relax()
if bugLev >= 1:
  print "test/hi/inputCif: === vasp ===\n%s\n=== end vasp ===" % (vasp,)

vasp.program = os.path.expanduser('~/pyladaExec.sh')   # default vasp program
vasp.has_nlep = True
setTestValidProgram( None)


vasp.isym       = 0
vasp.prec       = "accurate"
vasp.addgrid    = True
vasp.lmaxmix    = 4
vasp.lorbit     = 10

vasp.ismear     = 0
vasp.sigma      = 0.01
vasp.encut      = 340.0
vasp.nelm       = 1

vasp.lwave      = True
vasp.lmaxfockae = 4
vasp.nomega     = 64
##vasp.precfock   = False
vasp.encutgw    = 150.0
vasp.encutlf    = 150.0
vasp.lrpa       = False

vasp.isin       = 1
vasp.nbands     = 512
vasp.loptics    = True
vasp.lpead      = True

vasp.ldau       = True
##vasp.ldauprint  = 1
##vasp.ldautype   = 2
##vasp.ldul1      = 2 -1 -1
##vasp.lduu1      = -2.4 0.0 0.0
##vasp.lduj1      =  0.0 0.0 0.0
##vasp.lduo1      = 2 1 1



vasp.ediff        = 1e-3

vasp.npar       = 8
vasp.ncore      = 1
vasp.lplane     = True

vasp.set_symmetries = "off"
vasp.kpoints        = "\n0\nAuto\n20"

vasp.lvhar = False

vasp.lcharg = True


# Fast: for fast, low accuracy use:
# vasp.prec           = "low"
# vasp.ediff          = 1e-3
# vasp.encut          = 0.5
# vasp.kpoints        = "\n0\nGamma\n1 1 1\n0. 0. 0.\n"
# Or could try:
# vasp.kpoints        = "\n0\nAuto\n20\n"


# See vasp/functional.py:  elementName, fileName, max or min oxidation state

pseudoDir = "/nopt/nrel/ecom/cid/vasp.pseudopotentials.a/pseudos"

vasp.add_specie = "Ag", pseudoDir + "/Ag", U("dudarev", "d", 5.0), 1 
vasp.add_specie = "Al", pseudoDir + "/Al", None,  3
vasp.add_specie = "As", pseudoDir + "/As"
vasp.add_specie = "Au", pseudoDir + "/Au"
vasp.add_specie = "B",  pseudoDir + "/B"
vasp.add_specie = "Ba", pseudoDir + "/Ba", None, 2
vasp.add_specie = "Be", pseudoDir + "/Be"
vasp.add_specie = "Bi", pseudoDir + "/Bi"
vasp.add_specie = "Br", pseudoDir + "/Br"
vasp.add_specie = "C",  pseudoDir + "/C"
vasp.add_specie = "Ca", pseudoDir + "/Ca", None, 2
vasp.add_specie = "Cd", pseudoDir + "/Cd", U("dudarev", "d" , 6.0), 2
vasp.add_specie = "Cl", pseudoDir + "/Cl"
vasp.add_specie = "Co", pseudoDir + "/Co", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Cr", pseudoDir + "/Cr", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Cu", pseudoDir + "/Cu", U("dudarev", "d" , 5.0), 1
vasp.add_specie = "F",  pseudoDir + "/F"
vasp.add_specie = "Fe", pseudoDir + "/Fe", U("dudarev", "d", 3.0), 2
vasp.add_specie = "Ga", pseudoDir + "/Ga", None, 3
vasp.add_specie = "Ge", pseudoDir + "/Ge", None, 4
vasp.add_specie = "H",  pseudoDir + "/H",  None, 1
vasp.add_specie = "Hf", pseudoDir + "/Hf", U("dudarev", "d", 3.0), 2  # or U("dudarev", "f", 6.0), 2
vasp.add_specie = "Hg", pseudoDir + "/Hg"
vasp.add_specie = "I",  pseudoDir + "/I"
vasp.add_specie = "In", pseudoDir + "/In", None ,3
vasp.add_specie = "Ir", pseudoDir + "/Ir", U("dudarev", "d", 3.0), 3
vasp.add_specie = "K",  pseudoDir + "/K"
vasp.add_specie = "La", pseudoDir + "/La", U("dudarev", "d", 3.0), 2  # or U("dudarev", "f", 6.0), 2  
vasp.add_specie = "Li", pseudoDir + "/Li"
vasp.add_specie = "Lu", pseudoDir + "/Lu"
vasp.add_specie = "Mg", pseudoDir + "/Mg", None, 2
vasp.add_specie = "Mn", pseudoDir + "/Mn", U("dudarev", "d", 3.0 ), 2
vasp.add_specie = "Mo", pseudoDir + "/Mo", U("dudarev", "d", 3.0), 3
vasp.add_specie = "N",  pseudoDir + "/N"
vasp.add_specie = "Na", pseudoDir + "/Na"
vasp.add_specie = "Nb", pseudoDir + "/Nb", U("dudarev", "d", 3.0), 5
vasp.add_specie = "Ni", pseudoDir + "/Ni", U("dudarev", "d", 3.0), 2
vasp.add_specie = "O",  pseudoDir + "/O",  None, -2
vasp.add_specie = "P",  pseudoDir + "/P"
vasp.add_specie = "Pb", pseudoDir + "/Pb"
vasp.add_specie = "Pd", pseudoDir + "/Pd"
vasp.add_specie = "Pt", pseudoDir + "/Pt"
vasp.add_specie = "Rb", pseudoDir + "/Rb"
vasp.add_specie = "Rh", pseudoDir + "/Rh", U("dudarev", "d", 3.0), 3
vasp.add_specie = "Ru", pseudoDir + "/Ru"
vasp.add_specie = "S",  pseudoDir + "/S"
vasp.add_specie = "Sb", pseudoDir + "/Sb"
vasp.add_specie = "Sc", pseudoDir + "/Sc", U("dudarev", "d", 3.0),3
vasp.add_specie = "Se", pseudoDir + "/Se"
vasp.add_specie = "Si", pseudoDir + "/Si", None, 4
vasp.add_specie = "Sn", pseudoDir + "/Sn", None, 4
vasp.add_specie = "Sr", pseudoDir + "/Sr", None, 2
vasp.add_specie = "Ta", pseudoDir + "/Ta", U("dudarev", "d", 3.0),3
vasp.add_specie = "Te", pseudoDir + "/Te"
vasp.add_specie = "Ti", pseudoDir + "/Ti", U("dudarev", "d", 3.0), 4
vasp.add_specie = "V",  pseudoDir + "/V",  U("dudarev", "d", 3.0), 5
vasp.add_specie = "W",  pseudoDir + "/W",  U("dudarev", "d", 3.0), 5
vasp.add_specie = "Y",  pseudoDir + "/Y",  U("dudarev", "d", 3.0), 3
vasp.add_specie = "Zn", pseudoDir + "/Zn", U("dudarev", "d", 6.0), 2
vasp.add_specie = "Zr", pseudoDir + "/Zr", U("dudarev", "d", 3.0), 4


vasp.species["Co"].moment = [ 1.0, 4.0 ] 
vasp.species["Cr"].moment = [ 1.0, 4.0 ]
vasp.species["Cu"].moment = [ 1.0, 4.0 ] 
vasp.species["Fe"].moment = [ 1.0, 4.0 ]
vasp.species["Ir"].moment = [ 1.0, 4.0 ]
vasp.species["Mn"].moment = [ 1.0, 4.0 ]
vasp.species["Mo"].moment = [ 1.0, 4.0 ]
vasp.species["Ni"].moment = [ 1.0, 4.0 ]
vasp.species["Rh"].moment = [ 1.0, 4.0 ]
vasp.species["Sc"].moment = [ 1.0, 4.0 ]
vasp.species["Ti"].moment = [ 1.0, 4.0 ]
vasp.species["V" ].moment = [ 1.0, 4.0 ]


vasp.first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }
""" parameter to override during first relaxation step. """

vasp.relaxation = "relgw"
""" Degrees of freedom to relax. """

vasp.maxiter = 5
""" Maximum number of iterations before bailing out. """

vasp.keep_steps = True
""" Whether to keep or delete intermediate steps. """



