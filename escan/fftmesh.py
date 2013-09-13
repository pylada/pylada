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

""" Package to compute fft meshes from structural data. """
class SmallCells(object):
  """ Computes FFT mesh for small structures. 
  
      The FFT mesh is adjusted 
  """
  def __init__(self): 
    """ Initializes FFT mesh values. """
    super(SmallCells, self).__init__()

  def __call__(self, functional, structure, comm):
    """ Computes parameters of FFT mesh. """
    from operator import itemgetter
    from numpy import pi, sqrt
    from numpy.linalg import norm
    from pylada.physics import a0
    from quantities import angstrom

    para = structure.scale*2.0*sqrt(functional.cutoff)/pi/a0.rescale(angstrom).magnitude
    result = [int(norm(structure.cell[:,0]) * para + 5e-1), \
              int(norm(structure.cell[:,1]) * para + 5e-1), \
              int(norm(structure.cell[:,2]) * para + 5e-1)]

    assert result[0] * result[1] >= comm.size,\
           ValueError("Too many comms for a system this size.")

    if result[0] * result[1] % comm.size != 0:
      a = 0 if result[0] % comm.size  == 0 else comm.size - result[0] % comm.size
      b = 0 if result[1] % comm.size  == 0 else comm.size - result[1] % comm.size
      all_prods = [(i, j, i*j) for i in xrange(result[0], result[0]+a+1)\
                               for j in xrange(result[1], result[1]+b+1)\
                               if (i*j) % comm.size == 0]
      result[0], result[1], dummy = min(all_prods, key = itemgetter(2))

    return result, result, (0,0,0)

  def __repr__(self):
    """ Representation string of this object. """
    return "{0.__class__.__name__}()".format(self)


class Nanowire(object):
  """ Computes FFT mesh for Lijun's nanowires. """
  def ___init__(self):
    """ Initializes FFT mesh values. """
    super(SmallCells, self).__init__()

  def __call__(self, functional, structure, comm):
    """ Computes parameters of FFT mesh. """
    from operator import itemgetter
    from numpy import pi, sqrt
    from numpy.linalg import norm
    from pylada.physics import a0
    from quantities import angstrom

    para = structure.scale*2.0*sqrt(functional.cutoff)/pi/a0.rescale(angstrom).magnitude

    # first determines parameters of mesh.
    mesh = [int(norm(structure.cell[:,0]) * para + 5e-1), \
            int(norm(structure.cell[:,1]) * para + 5e-1), \
            int(norm(structure.cell[:,2]) * para + 5e-1)]

    multiple = comm.size if mesh[0] <= 500 else 3 * comm.size 
    if mesh[0] * mesh[1] % multiple != 0:
      a = 0 if mesh[0] % multiple  == 0 else multiple - mesh[0] % multiple
      b = 0 if mesh[1] % multiple  == 0 else multiple - mesh[1] % multiple
      all_prods = [(i, j, i*j) for i in xrange(mesh[0], mesh[0]+a+1)\
                               for j in xrange(mesh[1], mesh[1]+b+1)\
                               if (i*j) % multiple == 0]
      mesh[0], mesh[1], dummy = min(all_prods, key = itemgetter(2))

    # then determines parameters of smaller mesh.
    smesh = mesh if mesh[0] <= 500 else (mesh[0] // 3, mesh[1] // 3, mesh[2])
    if mesh[0] % smesh[0] != 0: smesh[0] += smesh[0] - mesh[0] % smesh[0]
    if mesh[1] % smesh[1] != 0: smesh[1] += smesh[1] - mesh[1] % smesh[1]

    # finally determines overlap mesh.
    omesh = (smesh[0]//2 if smesh[0] != mesh[0] else 0),\
            (smesh[1]//2 if smesh[1] != mesh[1] else 0),\
            (smesh[2]//2 if smesh[2] != mesh[2] else 0)

    return mesh, smesh, omesh

  def __repr__(self):
    """ Representation string of this object. """
    return "{0.__class__.__name__}()".format(self)
