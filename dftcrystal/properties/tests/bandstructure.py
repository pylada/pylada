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

def test():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import join
  from numpy import all, array
  from quantities import hartree
  from pylada.dftcrystal import Crystal, Functional, Shell
  from pylada.dftcrystal.properties import Properties
  from pylada import default_comm

  functional = Functional()
  functional.basis['Si'] = [
      Shell('s', a0=[16120.0, 0.001959],
                 a1=[2426.0, 0.01493], 
                 a2=[553.9, 0.07285],
                 a3=[156.3, 0.2461], 
                 a4=[50.07, 0.4859],
                 a5=[17.02, 0.325]),
      Shell('sp', a0=[292.7, -0.002781, 0.004438],
                  a1=[69.87, -0.03571, 0.03267],
                  a2=[22.34, -0.115, 0.1347],
                  a3=[8.15, 0.09356, 0.3287],
                  a4=[3.135, 0.603, 0.4496]), 
      Shell('sp', 4.0, a0=[1.22, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.55, 1.0, 1.0]),
      Shell('sp', 0.0, a0=[0.27, 1.0, 1.0]) ]

  functional.dft.pbe0 = True
  functional.fmixing = 30
  functional.shrink = 8, 16
  functional.levshift = 5, True
  functional.maxcycle = 600
  functional.dft.spin = True

  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  directory =  mkdtemp() 
  firstdir, seconddir = join(directory, '0'), join(directory, '1')
  try: 
     si = functional(crystal, outdir=firstdir, comm=default_comm)
     properties = Properties(si)
     properties.band = [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]
     result0 = properties(outdir=firstdir)
     assert result0.success 
     assert all(abs(result0.bandstructure.eigenvalues[0] - result0.bandstructure.eigenvalues[1]) < 1e-6)
     test = [ [-64.55637917, -64.55637522,  -3.64877203,  -3.64854929,
                -1.98183315,  -1.98183315,  -1.98183315,  -1.98097196,
                -1.98097196,  -1.98097196,   1.00052165,   1.49996607,
                 1.49996607,   1.49996607,   1.63212276,   1.63212276,
                 1.63212276,   1.66407077],
               [-64.55637917, -64.55637522,  -3.64877203,  -3.64854929,
                 -1.98183315,  -1.98183315,  -1.98183315,  -1.98097196,
                 -1.98097196,  -1.98097196,   1.00052165,   1.49996607,
                  1.49996607,   1.49996607,   1.63212276,   1.63212276,
                  1.63212276,   1.66407077] ] * hartree
     assert all(abs(result0.bandstructure.eigenvalues[0,0] - test[0]) < 1e-8)
     assert all(abs(result0.bandstructure.eigenvalues[0,-1] - test[1]) < 1e-8)
     test = array([[  6.66666667e-02,   0, 0],
                   [  1.33333333e-01,   0, 0],
                   [  2.00000000e-01,   0, 0],
                   [  2.66666667e-01,   0, 0],
                   [  3.33333333e-01,   0, 0],
                   [  4.00000000e-01,   0, 0],
                   [  4.66666667e-01,   0, 0],
                   [  5.33333333e-01,   0, 0],
                   [  6.00000000e-01,   0, 0],
                   [  6.66666667e-01,   0, 0],
                   [  7.33333333e-01,   0, 0],
                   [  8.00000000e-01,   0, 0],
                   [  8.66666667e-01,   0, 0],
                   [  9.33333333e-01,   0, 0],
                   [  1.00000000e+00,   0, 0],
                   [  1.06666667e+00,   0, 0],
                   [  9.58333333e-01,   4.16666667e-02,   0],
                   [  9.16666667e-01,   8.33333333e-02,   0],
                   [  8.75000000e-01,   1.25000000e-01,   0],
                   [  8.33333333e-01,   1.66666667e-01,   0],
                   [  7.91666667e-01,   2.08333333e-01,   0],
                   [  7.50000000e-01,   2.50000000e-01,   0],
                   [  7.08333333e-01,   2.91666667e-01,   0],
                   [  6.66666667e-01,   3.33333333e-01,   0],
                   [  6.25000000e-01,   3.75000000e-01,   0],
                   [  5.83333333e-01,   4.16666667e-01,   0],
                   [  5.41666667e-01,   4.58333333e-01,   0],
                   [  5.00000000e-01,   5.00000000e-01,   0],
                   [  4.58333333e-01,   5.41666667e-01,   0],
                   [  4.16666667e-01,   5.83333333e-01,   0],
                   [  3.75000000e-01,   6.25000000e-01,   0],
                   [  3.33333333e-01,   6.66666667e-01,   0],
                   [  2.91666667e-01,   7.08333333e-01,   0],
                   [  2.50000000e-01,   7.50000000e-01,   0],
                   [  2.08333333e-01,   7.91666667e-01,   0],
                   [  1.66666667e-01,   8.33333333e-01,   0],
                   [  1.25000000e-01,   8.75000000e-01,   0],
                   [  8.33333333e-02,   9.16666667e-01,   0],
                   [  4.16666667e-02,   9.58333333e-01,   0],
                   [ -1.11022302e-16,   1.00000000e+00,   0]])
     assert all(abs(result0.bandstructure.kpoints - test) < 1e-8)
     
     functional.dft.spin = False
     si = functional(crystal, outdir=seconddir, comm=default_comm)
     properties = Properties(si)
     properties.band = [0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0]
     result1 = properties(outdir=seconddir)
     assert result1.success 
     assert all(abs(result1.bandstructure.eigenvalues - result0.bandstructure.eigenvalues[0]) < 1e-6)
     assert all(abs(result1.bandstructure.kpoints - result0.bandstructure.kpoints) < 1e-6)
  finally: 
    if directory != '/tmp/test/': rmtree(directory)
if __name__ == '__main__':
  test()
