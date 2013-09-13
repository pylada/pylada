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

""" Tests that relaxed structure are correctly obtained. """
def compare_structures(a, b):
  """ compares two structures the long way. """
  from numpy import abs, all

  assert all(abs(a.cell - b.cell) < 1e-8)
  assert abs(a.scale - b.scale) < 1e-8 * a.scale.units
  assert set(a.__dict__.keys()) == set(b.__dict__.keys())
  for key, value in a.__dict__.iteritems():
    assert repr(value) == repr(getattr(b, key))
  assert len(a) == len(b)
  for atom0, atom1 in zip(a, b):
    assert all(abs(atom0.pos - atom1.pos) < 1e-8)
    assert atom0.type == atom1.type
    assert set(atom0.__dict__.keys()) == set(atom1.__dict__.keys())
    for key, value in atom0.__dict__.iteritems():
      assert repr(value) == repr(getattr(atom1, key))


def test_relax0(path):
  """ Tests structure with breaksym and relaxation with keepsymm. """
  from os.path import join
  from numpy import abs, array, all, dot
  from numpy.linalg import inv
  from quantities import e, hartree
  from pylada.dftcrystal import Extract
  from pylada.crystal import Structure

  extract = Extract(join(path, 'relax0.out'))
  # check some extraction stuff first
  assert extract.success
  assert extract._is_optgeom 
  assert extract.dimensionality == 3 
  assert abs(extract.total_energy + 5.786373101630E+02*hartree) < 1e-8
  
  # check-out that input structure is OK.
  assert len(extract.input_crystal) == 1
  assert abs(array(extract.input_crystal.params) - 5.43) < 1e-8
  assert extract.input_crystal.symmgroup == 227
  assert len(extract.input_crystal.atoms) == 1
  assert all(abs(extract.input_crystal.atoms[0].pos - 0.125) < 1e-8)
  assert extract.input_crystal.atoms[0].type == 'Si'
  assert extract.input_crystal[0].keyword.lower() == 'atomdisp'
  assert len(extract.input_crystal[0].atoms) == 1
  assert all(abs( extract.input_crystal[0].atoms[0].pos
                  - [0.01, 0, 0] ) < 1e-8)
  assert extract.input_crystal[0].atoms[0].type == 1
  assert not extract.input_crystal.is_bohr
  assert not extract.input_crystal.is_fractional
  assert extract.input_crystal.is_angstrom
  assert extract.input_crystal.is_breaksym

  # check that relaxation is as expected.
  assert extract.params.optgeom.enabled 
  assert extract.params.optgeom.keepsymm 

  # check that input structure is correctly extracted.
  structure = Structure( 0, 2.715, 2.715,
                         2.715, 0, 2.715,
                         2.715, 2.715, 0,
                         scale=1 )\
                .add_atom( 0.68875, 0.67875, 0.67875, 'Si',
                           asymmetric=True, group='SI', label=1)               \
                .add_atom( -0.67875, -0.67875, -0.67875, 'Si', 
                           asymmetric=True, group='SI', label=2)
  compare_structures(structure, extract.input_structure)
  compare_structures(structure, extract.input_crystal.eval())

  # check that output structure is correctly extracted.
  assert repr(extract.crystal[0]) == repr(extract.input_crystal[0])
  assert len(extract.crystal) == 4
  assert extract.crystal[1].keyword == 'keepsymm'
  assert extract.crystal[2].keyword == 'elastic'
  assert not extract.crystal[2].const_volume 
  assert extract.crystal[2].is_epsilon
  assert all(abs( extract.crystal[2].matrix 
                  - [ [-0.03907007718600364, 0, 0],
                      [0, -0.03906244198186859, -2.293447228964945e-06], 
                      [0, -2.293447228964945e-06, -0.03906244198186859] ]) < 1e-8)
  assert extract.crystal[3].keyword == 'atomdisp'
  assert len(extract.crystal[3].atoms) == 2
  assert all(abs(extract.crystal[3].atoms[0].pos - [-0.00478483, 0, 0]) < 1e-8)
  assert all(abs(extract.crystal[3].atoms[1].pos - [0.00478483, 0, 0]) < 1e-8)
  assert extract.crystal[3].atoms[0].type == 1
  assert extract.crystal[3].atoms[1].type == 2
  assert not extract.crystal.is_bohr
  assert not extract.crystal.is_fractional
  assert extract.crystal.is_angstrom
  assert not extract.crystal.is_breaksym

  structure = Structure( 
      0.000000000000E+00,  0.260892474044E+01,  0.260892474044E+01,
      0.260893924331E+01, -0.622670845305E-05,  0.260894547002E+01,
      0.260893924331E+01,  0.260894547002E+01, -0.622670845305E-05,
      scale=1 )                                                                \
      .add_atom( 6.570556585110E-01, 6.522348108275E-01, 6.522348108275E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e,
                  group='SI', label=1 )                                        \
      .add_atom( -6.474463592829E-01, -6.522348108275E-01, -6.522348108275E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e,
                 group='SI', label=2 )
  compare_structures(structure, extract.structure)
  for atom in structure:
    assert abs(atom.charge - 14*e) < 1e-8
    del atom.charge
  compare_structures(structure, extract.crystal.eval())

  posonly = extract._update_pos_only
  for atom, cmp in zip(posonly, structure): 
    pos = dot(dot(structure.cell, inv(posonly.cell)), atom.pos)
    assert all(abs(cmp.pos - pos) < 1e-8)

def test_relaxbhor(path):
  from os.path import join
  from numpy import abs, array, all, dot
  from numpy.linalg import inv
  from quantities import e, hartree
  from pylada.dftcrystal import Extract
  from pylada.crystal import Structure

  extract = Extract(join(path, 'relaxbhor.out'))
  # check some extraction stuff first
  assert extract.success
  assert extract._is_optgeom 
  assert extract.dimensionality == 3 
  assert abs(extract.total_energy + 5.786373123122E+02*hartree) < 1e-8
  
  # check-out that input structure is OK.
  assert len(extract.input_crystal) == 2
  assert abs(array(extract.input_crystal.params) - 5.43) < 1e-8
  assert extract.input_crystal.symmgroup == 227
  assert len(extract.input_crystal.atoms) == 1
  assert all(abs(extract.input_crystal.atoms[0].pos - 0.125) < 1e-8)
  assert extract.input_crystal.atoms[0].type == 'Si'
  assert extract.input_crystal[0].keyword.lower() == 'bohr'
  assert extract.input_crystal[1].keyword.lower() == 'atomdisp'
  assert len(extract.input_crystal[1].atoms) == 1
  assert all(abs(extract.input_crystal[1].atoms[0].pos - [0.01, 0, 0]) < 1e-8)
  assert extract.input_crystal[1].atoms[0].type == 1
  assert extract.input_crystal.is_bohr
  assert not extract.input_crystal.is_fractional
  assert not extract.input_crystal.is_angstrom
  assert extract.input_crystal.is_breaksym

  # check that relaxation is as expected.
  assert extract.params.optgeom.enabled 
  assert extract.params.optgeom.keepsymm 

  # check that input structure is correctly extracted.
  structure = Structure( 0, 2.715, 2.715,
                         2.715, 0, 2.715,
                         2.715, 2.715, 0,
                         scale=1 )\
                .add_atom( 0.684041772083, 0.67875, 0.67875, 'Si',
                           asymmetric=True, group='SI', label=1)               \
                .add_atom( -0.67875, -0.67875, -0.67875, 'Si', 
                           asymmetric=True, group='SI', label=2)
  compare_structures(structure, extract.input_structure)
  compare_structures(structure, extract.input_crystal.eval())

  # check that output structure is correctly extracted.
  assert len(extract.crystal) == 6
  for i in xrange(2):
    assert repr(extract.crystal[i]) == repr(extract.input_crystal[i])
  assert extract.crystal[2].keyword == 'keepsymm'
  assert extract.crystal[3].keyword == 'angstrom'
  assert extract.crystal[4].keyword == 'elastic'
  assert not extract.crystal[4].const_volume 
  assert extract.crystal[4].is_epsilon
  assert all(abs( extract.crystal[4].matrix 
                  - [ [-0.03900013662615098, 0, 0],
                      [0, -0.03900432708911383, -0.00010054972488253704], 
                      [0, -0.00010054972488253704, -0.03900432708911383] ]) < 1e-8)
  assert extract.crystal[5].keyword == 'atomdisp'
  assert len(extract.crystal[5].atoms) == 2
  assert all(abs(extract.crystal[5].atoms[0].pos - [-0.00252773, 0, 0]) < 1e-8)
  assert all(abs(extract.crystal[5].atoms[1].pos - [0.00252773, 0, 0]) < 1e-8)
  assert extract.crystal[5].atoms[0].type == 1
  assert extract.crystal[5].atoms[1].type == 2
  assert not extract.crystal.is_bohr
  assert not extract.crystal.is_fractional
  assert extract.crystal.is_angstrom
  assert not extract.crystal.is_breaksym

  structure = Structure( 
      0.000000000000E+00,  0.260911462906E+01,  0.260911462906E+01,
      0.260883025945E+01, -0.272992506112E-03,  0.260910325195E+01,
      0.260883025945E+01,  0.260910325195E+01, -0.272992506112E-03, 
      scale=1 )                                                                \
      .add_atom(  6.548363228522E-01,  6.522075648615E-01,  6.522075648615E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e,
                 group='SI', label=1)                                          \
      .add_atom( -6.497509306034E-01, -6.522075648615E-01, -6.522075648615E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e, 
                 group='SI', label=2)
  compare_structures(structure, extract.structure)

  for atom in structure:
    assert abs(atom.charge - 14*e) < 1e-8
    del atom.charge
  compare_structures(structure, extract.crystal.eval())
  posonly = extract._update_pos_only
  for atom, cmp in zip(posonly, structure): 
    pos = dot(dot(structure.cell, inv(posonly.cell)), atom.pos)
    assert all(abs(cmp.pos - pos) < 1e-8)

def test_relaxfrac(path):
  from os.path import join
  from numpy import abs, array, all, dot
  from numpy.linalg import inv
  from quantities import e, hartree
  from pylada.dftcrystal import Extract
  from pylada.crystal import Structure

  extract = Extract(join(path, 'relaxfrac.out'))
  # check some extraction stuff first
  assert extract.success
  assert extract._is_optgeom 
  assert extract.dimensionality == 3 
  assert abs(extract.total_energy + 5.786373124870E+02*hartree) < 1e-8

  # check-out that input structure is OK.
  assert len(extract.input_crystal) == 3
  assert abs(array(extract.input_crystal.params) - 5.43) < 1e-8
  assert extract.input_crystal.symmgroup == 227
  assert len(extract.input_crystal.atoms) == 1
  assert all(abs(extract.input_crystal.atoms[0].pos - 0.125) < 1e-8)
  assert extract.input_crystal.atoms[0].type == 'Si'
  assert extract.input_crystal[0].keyword.lower() == 'fraction'
  assert extract.input_crystal[1].keyword.lower() == 'keepsymm'
  assert extract.input_crystal[2].keyword.lower() == 'atomdisp'
  assert len(extract.input_crystal[2].atoms) == 1
  assert all(abs(extract.input_crystal[2].atoms[0].pos - [0.01, 0, 0]) < 1e-8)
  assert extract.input_crystal[2].atoms[0].type == 1
  assert not extract.input_crystal.is_bohr
  assert extract.input_crystal.is_fractional
  assert not extract.input_crystal.is_angstrom
  assert not extract.input_crystal.is_breaksym

  # check that relaxation is as expected.
  assert extract.params.optgeom.enabled 
  assert extract.params.optgeom.keepsymm 

  # check that input structure is correctly extracted.
  structure = Structure( 0, 2.715, 2.715,
                         2.715, 0, 2.715,
                         2.715, 2.715, 0,
                         scale=1 )\
                .add_atom( 0.67875, 0.7059, 0.7059, 'Si',
                           asymmetric=True, group='SI', label=1)               \
                .add_atom( -0.7059, -0.67875, -0.7059, 'Si', 
                           asymmetric=False, group='SI', label=2)
  compare_structures(structure, extract.input_structure)
  compare_structures(structure, extract.input_crystal.eval())
  
  # check that output structure is correctly extracted.
  for i in xrange(3):
    assert repr(extract.crystal[i]) == repr(extract.input_crystal[i])
  assert len(extract.crystal) == 6
  assert extract.crystal[3].keyword == 'angstrom'
  assert extract.crystal[4].keyword == 'elastic'
  assert not extract.crystal[4].const_volume 
  assert extract.crystal[4].is_epsilon
  assert all(abs( extract.crystal[4].matrix 
       - [ [-0.03909399840733285, 1.9826783023479955e-05, -7.783463753630038e-06],
           [1.9826783023479955e-05, -0.03909399840733285, -7.783463753630038e-06],
           [-5.8938764164088395e-06, -5.8938764164088395e-06, -0.039119955865756784] ]) < 1e-8)
  assert extract.crystal[5].keyword == 'atomdisp'
  assert len(extract.crystal[5].atoms) == 1
  assert all(abs( extract.crystal[5].atoms[0].pos
                  - [-0.013094068944795376, -0.01309406894477713, -0.02598436677348225] ) < 1e-8)
  assert extract.crystal[5].atoms[0].type == 1
  assert not extract.crystal.is_bohr
  assert not extract.crystal.is_fractional
  assert extract.crystal.is_angstrom
  assert not extract.crystal.is_breaksym

  structure = Structure( 
      0.326976118176E-04, 0.260883866222E+01,  0.260891362404E+01,
      0.260883866222E+01, 0.326976118176E-04,  0.260891362404E+01,
      0.260877331795E+01, 0.260877331795E+01, -0.320037489411E-04,
      scale=1 )                                                                \
      .add_atom(  6.391293810153E-01,  6.652174406614E-01,  6.522926954249E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e,
                 group='SI', label=1)                                          \
      .add_atom( -6.652174406614E-01, -6.391293810153E-01, -6.522926954249E-01,
                 'Si', asymmetric=False, charge=array(14.0) * e, 
                 group='SI', label=2)
  compare_structures(structure, extract.structure)

  for atom in structure:
    assert abs(atom.charge - 14*e) < 1e-8
    del atom.charge
  compare_structures(structure, extract.crystal.eval())
  posonly = extract._update_pos_only
  for atom, cmp in zip(posonly, structure): 
    pos = dot(dot(structure.cell, inv(posonly.cell)), atom.pos)
    assert all(abs(cmp.pos - pos) < 1e-8)

def test_relaxfrac2(path):
  from os.path import join
  from numpy import abs, array, all, dot
  from numpy.linalg import inv
  from quantities import e, hartree
  from pylada.dftcrystal import Extract
  from pylada.crystal import Structure

  extract = Extract(join(path, 'relaxfrac2.out'))
  # check some extraction stuff first
  assert extract.success
  assert extract._is_optgeom 
  assert extract.dimensionality == 3 
  assert abs(extract.total_energy + 5.786359136346E+02*hartree) < 1e-8

  # check-out that input structure is OK.
  assert len(extract.input_crystal) == 6
  assert abs(array(extract.input_crystal.params) - 5.43) < 1e-8
  assert extract.input_crystal.symmgroup == 227
  assert len(extract.input_crystal.atoms) == 1
  assert all(abs(extract.input_crystal.atoms[0].pos - 0.125) < 1e-8)
  assert extract.input_crystal.atoms[0].type == 'Si'
  assert extract.input_crystal[0].keyword.lower() == 'fraction'
  assert extract.input_crystal[1].keyword.lower() == 'keepsymm'
  assert extract.input_crystal[2].keyword.lower() == 'atomdisp'
  assert len(extract.input_crystal[2].atoms) == 1
  assert all(abs(extract.input_crystal[2].atoms[0].pos - [0.01, 0, 0]) < 1e-8)
  assert extract.input_crystal[2].atoms[0].type == 1
  assert extract.input_crystal[3].keyword.lower() == 'angstrom' 
  assert extract.input_crystal[4].keyword == 'elastic'
  assert not extract.input_crystal[4].const_volume 
  assert extract.input_crystal[4].is_epsilon
  assert all(abs( extract.input_crystal[4].matrix 
       - [ [-0.0390939984073, 1.98267830235e-05, -7.78346375363e-06],
           [1.98267830235e-05, -0.0390939984073, -7.78346375363e-06],
           [-5.89387641641e-06, -5.89387641641e-06, -0.0391199558658] ]) < 1e-8)
  assert extract.input_crystal[5].keyword == 'atomdisp'
  assert len(extract.input_crystal[5].atoms) == 1
  assert all(abs( extract.input_crystal[5].atoms[0].pos
                  - [-0.013094068944795376, -0.01309406894477713, -0.02598436677348225] ) < 1e-8)
  assert extract.input_crystal[5].atoms[0].type == 1
  assert not extract.input_crystal.is_bohr
  assert not extract.input_crystal.is_fractional
  assert extract.input_crystal.is_angstrom
  assert not extract.input_crystal.is_breaksym

  # check that relaxation is as expected.
  assert extract.params.optgeom.enabled 
  assert extract.params.optgeom.keepsymm 

  # check that input structure is correctly extracted.
  structure = Structure(
        0.326976118177E-04, 0.260883866222E+01,  0.260891362404E+01,
        0.260883866222E+01, 0.326976118177E-04,  0.260891362404E+01,
        0.260877331795E+01, 0.260877331795E+01, -0.320037489411E-04,
        scale=1 )\
        .add_atom( 6.391293810153E-01,  6.652174406614E-01,  6.522926954249E-01,
                   'Si', asymmetric=True, group='SI', label=1)\
        .add_atom(-6.652174406614E-01, -6.391293810153E-01, -6.522926954249E-01,
                   'Si', asymmetric=False, group='SI', label=2)
  compare_structures(structure, extract.input_structure)
  compare_structures(structure, extract.input_crystal.eval())
  
  # check that output structure is correctly extracted.
  assert len(extract.crystal) == 8
  for i in xrange(6):
    assert repr(extract.crystal[i]) == repr(extract.input_crystal[i])
  assert extract.crystal[6].keyword == 'elastic'
  assert extract.crystal[6].is_epsilon
  assert all(abs( extract.crystal[6].matrix 
       - [ [0.0014910603539926015, -3.261366539075716e-05, 9.243596223229567e-06],
           [-3.2613665390868185e-05, 0.0014910603539930456, 9.243596223229567e-06],
           [9.24445028060031e-06, 9.244450281044398e-06, 0.0015205063395626883] ]) < 1e-8)
  assert extract.crystal[7].keyword == 'atomdisp'
  assert len(extract.crystal[7].atoms) == 1
  assert all(abs( extract.crystal[7].atoms[0].pos
                  - [4.5230400922175564e-05, 4.523040092304581e-05, -7.99993002858215e-05] ) < 1e-8)
  assert extract.crystal[7].atoms[0].type == 1
  assert not extract.crystal.is_bohr
  assert not extract.crystal.is_fractional
  assert extract.crystal.is_angstrom
  assert not extract.crystal.is_breaksym

  structure = Structure( 
      -0.282229780687E-04,  0.261275271150E+01, 0.261271858518E+01,
       0.261275271150E+01, -0.282229780687E-04, 0.261271858518E+01,
       0.261276409190E+01,  0.261276409190E+01, 0.161835337242E-04, 
      scale=1 )                                                                \
      .add_atom(  6.401119262486E-01,  6.662397355934E-01,  6.532165692726E-01,
                 'Si', asymmetric=True, charge=array(14.0) * e,
                 group='SI', label=1)                                          \
      .add_atom( -6.662397355934E-01, -6.401119262486E-01, -6.532165692726E-01,
                 'Si', asymmetric=False, charge=array(14.0) * e, 
                 group='SI', label=2)
  compare_structures(structure, extract.structure)

  for atom in structure:
    assert abs(atom.charge - 14*e) < 1e-8
    del atom.charge
  compare_structures(structure, extract.crystal.eval())
  posonly = extract._update_pos_only
  for atom, cmp in zip(posonly, structure): 
    pos = dot(dot(structure.cell, inv(posonly.cell)), atom.pos)
    assert all(abs(cmp.pos - pos) < 1e-8)

if __name__ == '__main__':
  from sys import argv
  from os.path import join, dirname

  test_relax0(join(dirname(argv[0]), 'data'))
  test_relaxbhor(join(dirname(argv[0]), 'data'))
  test_relaxfrac(join(dirname(argv[0]), 'data'))
  test_relaxfrac2(join(dirname(argv[0]), 'data'))
