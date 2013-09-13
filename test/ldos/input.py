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

# sets up input stuff
vff = Vff()
vff.lattice.set_types = ("Si", "Ge"), ("Si","Ge")
vff.lattice.scale = 5.45
vff.add_bond = "Si", "Si", (2.35905923, 45.375785351, -150.363195485, 519.230157133)
vff.add_bond = "Ge", "Ge", (2.44167023, 37.816422666, -124.830189020, 250.179861133)
vff.add_bond = "Si", "Ge", (2.39642606, 41.875535703, -145.994430120, 366.853558523)
vff.add_angle = "Si", "Si", "Si", ("tet", 0, 12.579328566, 44.930324684, 359.422663897)
vff.add_angle = "Ge", "Ge", "Ge", ("tet", 0, 10.383093493, 55.677487481)
vff.add_angle = "Si", "Ge", "Si", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Ge", "Si", "Ge", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Si", "Si", "Ge", ("tet", 0, 12.245971457, 32.561864518, 179.706525419)
vff.add_angle = "Si", "Ge", "Ge", ("tet", 0, 11.147853920, 37.930591733)
vff.minimizer.verbose = True
vff.minimizer.type = "frprmn"
vff.minimizer.itermax = 4000
vff.minimizer.tolerance = 1e-5
vff.minimizer.uncertainties = 1e-3

escan = Escan()
escan.eref                  = None
escan.cutoff                = 8.2
escan.smooth                = 1.0
escan.kinetic_scaling       = 1.0
escan.nbstates              = 10
escan.itermax               = 60
escan.nlines                = 200
escan.tolerance             = 1e-12
escan.rspace_cutoff         = 3.75
escan.fft_mesh              = fftmesh.SmallCells()
escan.do_genpot             = True
escan.do_relax              = True
escan.input_wavefunctions   = None
escan.kpoint                = array([0,0,0])
escan.dnc_mesh              = None
escan.overlap_mesh          = None
escan.potential             = soH
escan.add_potential         = "test_input/pseudos/vq.Si", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.0495
escan.add_potential         = "test_input/pseudos/vq.Ge", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.2774
escan.add_potential         = "test_input/pseudos/vq.SiGe.Si", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.0495
escan.add_potential         = "test_input/pseudos/vq.SiGe.Ge", "test_input/pseudos/vwr.pso",\
                              0,0,0,0.2774
escan.maskr                 = "test_input/pseudos/maskr"
escan.vff = vff
