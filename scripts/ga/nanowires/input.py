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

from pylada.escan import fftmesh
# sets up input stuff
vff = Vff()
vff.lattice.set_types = ("Si", "Ge", "Hg"), ("Si","Ge","Hg")
vff.add_bond = "Si", "Si", (2.35905923, 45.375785351, -150.363195485, 519.230157133)
vff.add_bond = "Ge", "Ge", (2.44167023, 37.816422666, -124.830189020, 250.179861133)
vff.add_bond = "Si", "Ge", (2.39642606, 41.875535703, -145.994430120, 366.853558523)
vff.add_bond = "Hg", "Hg", (2.35905923,  7.000000000)
vff.add_bond = "Si", "Hg", (2.35905923,  7.000000000)
vff.add_bond = "Ge", "Hg", (2.44167023,  7.000000000)
vff.add_angle = "Si", "Si", "Si", ("tet", 0, 12.579328566, 44.930324684, 359.422663897)
vff.add_angle = "Ge", "Ge", "Ge", ("tet", 0, 10.383093493, 55.677487481)
vff.add_angle = "Si", "Ge", "Si", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Ge", "Si", "Ge", ("tet", 0, 11.902809727, 20.193404352)
vff.add_angle = "Si", "Si", "Ge", ("tet", 0, 12.245971457, 32.561864518, 179.706525419)
vff.add_angle = "Si", "Ge", "Ge", ("tet", 0, 11.147853920, 37.930591733)
vff.add_angle = "Hg", "Hg", "Hg", ("tet", 0,  3.000000000)
vff.add_angle = "Si", "Hg", "Si", ("tet", 0,  3.000000000)
vff.add_angle = "Ge", "Hg", "Ge", ("tet", 0,  3.000000000)
vff.add_angle = "Si", "Hg", "Hg", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Hg", "Ge", ("tet", 0,  3.000000000)
vff.add_angle = "Ge", "Hg", "Si", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Si", "Hg", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Si", "Si", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Si", "Ge", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Ge", "Hg", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Ge", "Ge", ("tet", 0,  3.000000000)
vff.add_angle = "Hg", "Ge", "Si", ("tet", 0,  3.000000000)
vff.minimizer.verbose = True
vff.minimizer.type = "gsl_bfgs2"
vff.minimizer.itermax = 6000
vff.minimizer.tolerance = 1e-2
vff.minimizer.uncertainties = 1e-3

escan = BandGap()
escan.vff = vff
escan.cutoff                = 8.2
escan.smooth                = 1.0
escan.kinetic_scaling       = 1.0
escan.nbstates              = 10
escan.itermax               = 60
escan.nlines                = 200
escan.tolerance             = 1e-12
escan.rspace_cutoff         = 3.75
escan.fft_mesh              = fftmesh.Nanowire()
escan.do_genpot             = True
escan.do_relax              = True
escan.input_wavefunctions   = None
escan.kpoint                = array([0,0,0])
escan.dnc_mesh              = None
escan.overlap_mesh          = None
escan.potential             = soH
escan.maskr                 = "$HOME/nanowires/pseudos/maskr"
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Ge32032",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.2774
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Ge32014",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.2774
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Ge32080",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.2774
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Si14014",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0475
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Si14032",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0475
escan.add_potential         = "$HOME/nanowires/pseudos/vq.Si14080",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0475
escan.add_potential         = "$HOME/nanowires/pseudos/vq.psdH80080",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0495
escan.add_potential         = "$HOME/nanowires/pseudos/vq.psdH80014",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0495
escan.add_potential         = "$HOME/nanowires/pseudos/vq.psdH80032",\
                              "$HOME/nanowires/pseudos/vwr.pso",\
                              0.0,0.0,0.0,0.0495
escan.references            = None
escan.nbstates              = (2, 4)
degeneracy                  = 1e-3 # in eV


growth_direction = (0,0,1)
""" Growth direction. """
thickness = 0.35
""" Shell thickness. """
core_radius = 2
""" Radius of the nanowire core. """
core_types = ['Si', 'Ge']
""" Type of the nanowire core. """
passivant = 'Hg'
""" Type of the passivant. """

nb_trials = 4
""" Number of different trial GA runs. """

Evaluator = DipoleEvaluator
""" Kind of objective funtion we will use. """

ranges = [ (0, 8) ]
""" nanowire sizes. """

crossover_rate = 0.8
""" Rate of crossover vs other operations. """
swap_rate = 0.1
""" Rate of swap-type mutations over other operations. """
growth_rate = 0.1
""" Rate of growth/shrink-type mutations over other operations. """

offspring_rate = 0.2
""" Rate of new offspring creation per generation. """
population_size = 10
""" Size of the population. """

max_generations = -1
""" Maximum number of generation. """

pools = 5
""" Number of pools of processor with which to perform calculations. """

rootworkdir = "$SCRATCH"
""" Root of the directory where calculations are performed. 
   
    If None, then calculations and ouput are performed in the same directory tree.
    If something like "$SCRATCH", then calculations are perfomed in a directory
    where the root is defined by the environment variable SCRATCH, whereas
    output would be on a similar directory tree rooted in, most likely, the
    user's home directory.
"""
