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
  from numpy import all, abs, dot, identity
  from pylada.dftcrystal import Crystal, External
  from pylada.crystal import which_site

  orig = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  b = External(copy=orig.eval())
  assert all(abs(orig.eval().cell - b.initial.cell) < 1e-8)
  assert len(b.initial) == 2
  assert repr(b) == repr( eval(repr(b), {'External': External}) )
  assert len(b.print_input().splitlines()) == 2
  assert b.print_input().splitlines()[0] == 'EXTERNAL'
  assert b.print_input().splitlines()[1] == 'END EXTERNAL'

  # check that symmetries as extracted are ok.
  for atom in b.initial:
    for op in orig.symmetry_operators:
      assert all(abs(identity(3, dtype='float64') - dot(op[:3].T, op[:3])) < 1e-8)
      assert which_site(dot(op[:3], atom.pos) + op[3], b.initial) != -1

  # Check that running structure through crystal works
  c = b.eval()
  assert all(abs(c.cell - b.cell) < 1e-8)
  assert len(c) == 2
  assert all(abs(abs(c[0].pos) - [0.67875]*3) < 1e-8)
  assert all(abs(abs(c[1].pos) - [0.67875]*3) < 1e-8)
  assert all(abs(c[0].pos + c[1].pos) < 1e-8)

def test_readwrite():
  """ Test that writing to and reading from file yield same result. """
  from tempfile import mkdtemp
  from shutil import rmtree
  from pylada.dftcrystal import External, Supercell, InsertAtoms
  from pylada.crystal import read, write
  from pylada.misc import Changedir

  structure = External( 2.97938, 0, 0,
                        0, 6.56068, 0,
                        0, 0, 34.6664,
                        scale=1 )\
        .add_atom(0, -1.64017, 6.33322, 'O', asymmetric=True, group='O', label=1)\
        .add_atom(0, 1.64017, 4.7745, 'Ti', asymmetric=True, group='TI', label=2)\
        .add_atom(1.48969, -1.64017, 5.28628, 'Ti', asymmetric=True, group='TI', label=3)\
        .add_atom(1.48969, 2.872, 5.15099, 'O', asymmetric=True, group='O', label=4)\
        .add_atom(1.48969, 0.408334, 5.15099, 'O', asymmetric=False, group='O', label=5)\
        .add_atom(0, -1.64017, 3.768, 'O', asymmetric=True, group='O', label=6)\
        .add_atom(0, 1.64017, 2.98235, 'O', asymmetric=True, group='O', label=7)\
        .add_atom(0, -1.64017, 1.96284, 'Ti', asymmetric=True, group='TI', label=8)\
        .add_atom(1.48969, 1.64017, 1.51848, 'Ti', asymmetric=True, group='TI', label=9)\
        .add_atom(1.48969, -0.36011, 1.71256, 'O', asymmetric=True, group='O', label=10)\
        .add_atom(1.48969, -2.92023, 1.71256, 'O', asymmetric=False, group='O', label=11)\
        .add_atom(0, 1.64017, 0.407607, 'O', asymmetric=True, group='O', label=12)\
        .add_atom(0, -1.64017, -0.407607, 'O', asymmetric=False, group='O', label=13)\
        .add_atom(0, 1.64017, -1.96284, 'Ti', asymmetric=False, group='TI', label=14)\
        .add_atom(1.48969, -1.64017, -1.51848, 'Ti', asymmetric=False, group='TI', label=15)\
        .add_atom(1.48969, 2.92023, -1.71256, 'O', asymmetric=False, group='O', label=16)\
        .add_atom(1.48969, 0.36011, -1.71256, 'O', asymmetric=False, group='O', label=17)\
        .add_atom(0, -1.64017, -2.98235, 'O', asymmetric=False, group='O', label=18)\
        .add_atom(0, 1.64017, -3.768, 'O', asymmetric=False, group='O', label=19)\
        .add_atom(0, -1.64017, -4.7745, 'Ti', asymmetric=False, group='TI', label=20)\
        .add_atom(1.48969, 1.64017, -5.28628, 'Ti', asymmetric=False, group='TI', label=21)\
        .add_atom(1.48969, -0.408334, -5.15099, 'O', asymmetric=False, group='O', label=22)\
        .add_atom(1.48969, -2.872, -5.15099, 'O', asymmetric=False, group='O', label=23)\
        .add_atom(0, 1.64017, -6.33322, 'O', asymmetric=False, group='O', label=24)\
        .append(Supercell([[2, 0, 0], [0, 2, 0], [0, 0, 1]]))\
        .append(InsertAtoms()\
                  .add_atom(0, -1.64017, 7.33322, 'H'))

  directory = mkdtemp()
  try:
    with Changedir(directory) as cwd:
      write.crystal(structure.initial, 'fort.34')
      initial = read.crystal('fort.34')

      assert write.crystal(structure.initial, None) == write.crystal(initial, None)
  finally:
    try: rmtree(directory)
    except: pass

if __name__ == '__main__':
  test()
  test_readwrite()
