def test_smith():
  from numpy import dot, all
  from numpy.random import randint
  from numpy.linalg import det
  from pylada.math import smith_normal_form
  from pylada.error import input
  
  for i in xrange(50):
    cell = randint(-5, 5, size=(3,3))
    while abs(det(cell)) < 1e-2: cell = randint(-5, 5, size=(3,3))
    s, l, r = smith_normal_form(cell)
    assert all(dot(dot(l, cell), r) == s)

  try: smith_normal_form([[0, 0, 0], [1, 2, 0], [3, 4, 5]])
  except input: pass
  else: raise Exception()

def test_capi():
  from _gruber import testme
  assert testme()

if __name__ == '__main__':
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  test_capi()
  test_smith()
