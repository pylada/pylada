def test_crystal():
  from numpy import array
  from lada.dftcrystal.crystal import Crystal
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O')
  assert structure.print_input() is not None
  assert len(structure.print_input().split('\n')) == 9
  assert all(abs(array(structure.print_input().split()[1:-2], dtype='float64') \
             - [0, 0, 0, 136, 4.63909875, 2.97938395, 2, 22, 0, 0, 0, 8, 0.306153, 0.306153, 0]) < 1e-8)
  assert structure.print_input().split()[0] == 'CRYSTAL'
  assert structure.print_input().split()[-2] == 'END'
  assert structure.print_input().split()[-1] == 'CRYSTAL'

def test_supercell():
  from numpy import array
  from lada.dftcrystal import Supercell, Crystal
  from lada.error import ValueError
  structure = Crystal(136, 4.63909875, 2.97938395, \
                      ifhr=0, \
                      shift=0)\
                  .add_atom(0, 0, 0, 'Ti')\
                  .add_atom(0.306153, 0.306153, 0, 'O') \
                  .append(Supercell([[2,0,0],[0,1,0],[0,0,1]]))
  assert len(structure.eval()) == 2 * 6
  structure[-1].matrix[2,2] = 2
  assert len(structure.eval()) == 2 * 2 * 6

  # test errors
  try: structure[-1].matrix = array( [[0]*3]*3 )
  except ValueError: pass
  else: raise Exception()
  structure[-1].matrix = array([[1,0],[0,1]])
  try: structure.print_input()
  except ValueError: pass
  else: raise Exception()


def test_elastic():
  """ Test elastic keyword. """
  from numpy import array, all, abs, identity, dot
  from lada.dftcrystal import Crystal, Elastic

  epsilon = array([[0, 0, 0], [0, 0, 0], [0, 0, 0.1]])
  crystal = Crystal(227, 5.43).add_atom(0.125, 0.125, 0.125, 'Si')
  structure = crystal.eval()
  crystal.append(Elastic([[0, 0, 0], [0, 0, 0], [0, 0, 0.1]]))
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  assert abs(crystal.eval()[0].pos[1] - crystal.eval()[0].pos[0]) < 1e-8
  assert abs(crystal.eval()[1].pos[2] / crystal.eval()[1].pos[1] - 1.1) < 1e-8
  assert abs(crystal.eval()[1].pos[1] - crystal.eval()[1].pos[0]) < 1e-8
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)

  z = dot(epsilon, structure.cell)
  crystal[-1].matrix = z
  crystal[-1].is_epsilon = False
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)

  epsilon = array([[0, 0, 0], [0, 0, 0], [0, 0.1, 0]])
  crystal[-1].matrix = epsilon
  crystal[-1].is_epsilon = True
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)

  z = dot(epsilon, structure.cell)
  crystal[-1].matrix = z
  crystal[-1].is_epsilon = False
  assert all(abs(dot(identity(3) + epsilon, structure.cell) - crystal.eval().cell) < 1e-8)
  assert all(abs(crystal.eval()[0].pos - dot(identity(3) + epsilon, structure[0].pos)) < 1e-8)
  assert all(abs(crystal.eval()[1].pos - dot(identity(3) + epsilon, structure[1].pos)) < 1e-8)
  a = Elastic()
  a.raw = crystal[-1].raw
  assert a.is_epsilon == crystal[-1].is_epsilon
  assert all(abs(a.matrix - crystal[-1].matrix) < 1e-8)



if __name__ == '__main__':
  test_crystal()
  test_supercell()
  test_elastic()
