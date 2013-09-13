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

""" Input script to oscillator strength GA. """
# vff parameters
vff = Vff()
vff.lattice.set_types = ("Si", "Ge"), ("Si","Ge")
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
vff.minimizer.type = "gsl_bfgs2"
vff.minimizer.itermax = 4000
vff.minimizer.tolerance = 1e-5
vff.minimizer.uncertainties = 1e-3
vff.direction = FreezeCell.a1 | FreezeCell.a2

# escan parameters
escan = Escan()
escan.cutoff                = 8.2
escan.smooth                = 1.0
escan.kinetic_scaling       = 1.0
escan.nbstates              = 10
escan.itermax               = 20
escan.nlines                = 200
escan.tolerance             = 1e-6
escan.rspace_cutoff         = 3.75
escan.do_genpot             = True
escan.do_relax              = True
escan.input_wavefunctions   = None
escan.kpoint                = array([0,0,0])
escan.dnc_mesh              = None
escan.overlap_mesh          = None
escan.potential             = soH
escan.add_potential         = "~/SiGe/pseudos/vq.Si",      "~/SiGe/pseudos/vwr.pso", 0,0,0,0.0495
escan.add_potential         = "~/SiGe/pseudos/vq.Ge",      "~/SiGe/pseudos/vwr.pso", 0,0,0,0.2774
escan.add_potential         = "~/SiGe/pseudos/vq.SiGe.Si", "~/SiGe/pseudos/vwr.pso", 0,0,0,0.0495
escan.add_potential         = "~/SiGe/pseudos/vq.SiGe.Ge", "~/SiGe/pseudos/vwr.pso", 0,0,0,0.2774
escan.maskr                 = "~/SiGe/pseudos/maskr"
escan.vff = vff

# GA parameters.
population_size = 100
""" Size of the GA population """
max_generations = -1
""" Maximum number of generations """
offspring_rate  = 0.1
""" Rate at which offspring are created. """
crossover_vs_mutation_rate = 0.8
""" Rate of crossover operations over mutation operations. """
mean_conc = 0.2
""" Mean concentration of randomly created individual. 

    If None, use uniform distribution, rather than guaussian distribution.
"""
stddev_conc = 0.4
""" Standard deviation of the concentration of randomly created individual. """
pools = 5
""" Number of pools of processes.

    During the evaluation of new individuals, the processes are splitted
    amongst this many pools. Each individual is computed by a single pool.
    If there are restrictions (as for pescan) as too many procs can be used by
    each calculation, they  are for the user to figure out.
"""

# Configuration space.
def supercell(n): 
  """ Supercell within which to perform calculations. """
  result = array( [ [1, 0,    0],
                    [0, 1, -0.5],
                    [0, 1,  0.5] ], dtype="float64" )
  result[:, 0] *= float(n)
  return result

def fftmesh(cell):
  """ FFT real-space mesh used by escan. """
  result = int(norm(cell[:,0]) * 20e0 + 5e-1), \
           int(norm(cell[:,1]) * 20e0 + 5e-1), \
           int(norm(cell[:,2]) * 20e0 + 5e-1) 
  # This line makes sure mesh is even. Not sure why...
  result = result[0] + result[0] % 2, \
           result[1] + result[1] % 2, \
           result[2] + result[2] % 2
  return result

scales = [5.45]
""" Substrate lattice parameter. """
trials = range(2)
""" Number of trials per configuration space. """
periods = [20]
""" Size of the period in the configuration space. """

references = (-0.45, 0.10)
""" References for bandgap calculations. """
kwargs = { 'overlap_factor': 6e0, 'nbstates': 6 }
""" Extra arguments passed on to functional. """
