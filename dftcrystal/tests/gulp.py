def test():
  """ Tests writing out gulp files. """
  from numpy import array, all, abs
  from pylada.dftcrystal import Crystal
  from pylada.crystal.write import gulp

  a = Crystal(136, 4.63909875, 2.97938395, \
              ifhr=0, \
              shift=0)\
             .add_atom(0, 0, 0, 'Ti')\
             .add_atom(0.306153, 0.306153, 0, 'O')

  string = gulp(a).splitlines()
  assert string[0] == 'vectors'
  assert all(abs(array(string[1].split(), dtype='float64') - [4.63909875, 0, 0]) < 1e-5) 
  assert all(abs(array(string[2].split(), dtype='float64') - [0, 4.63909875, 0]) < 1e-5) 
  assert all(abs(array(string[3].split(), dtype='float64') - [0, 0, 2.97938395]) < 1e-5) 
  assert string[4] == 'spacegroup'
  assert string[5] == '136'
  assert string[6].rstrip().lstrip() == 'cartesian'
  assert string[7].split()[:2] == ['Ti', 'core']
  assert all(abs(array(string[7].split()[2:], dtype='float64')) < 1e-5)
  assert string[8].split()[:2] == ['O', 'core']
  assert all(abs(array(string[8].split()[2:], dtype='float64') - [1.420274, 1.420274, 0]) < 1e-5)

if __name__ == '__main__': 
  test()





