def test_isinteger():
  from numpy import array
  from lada.enum.cppwrappers import is_integer
  assert not is_integer(array([0.5, 0.5]))
  assert not is_integer(array([1.0, 0.5]))
  assert not is_integer(array(0.1))
  assert is_integer(array([1.0, -1.0])) 
  assert is_integer(array(-5))
  assert is_integer(array([ [1.0, -1.0], [5.0, 1111.0] ]))

def test_ndimiterator():
  from numpy import all 
  from itertools import product
  from lada.enum.cppwrappers import NDimIterator
  from lada.error import TypeError, ValueError

  i = 0
  for u in NDimIterator(3,3,4): i += 1; continue
  assert i == 3*3*4
  assert all(u == 1) # went round
  for u in NDimIterator(3,3,4): continue
  assert all(u == 1) # went round
  iterator = NDimIterator(5)
  for v in xrange(1, 6):
    u = iterator.next()
    assert  u == v
  assert u == 5 # didn't go round cos xrange finished first.
  iterator = NDimIterator(2, 3)
  for v in product(xrange(1, 3), xrange(1, 4)):
    u = iterator.next()
    assert all(u == v)
  assert all(u == [2, 3])
  iterator = NDimIterator(2, 3)
  a = iterator.next()
  assert all(a == [1, 1])
  b = iterator.next()
  assert all(b == [1, 2])
  assert all(a == b)
  assert a is b
  for u, v in zip( NDimIterator(5, 5, 5), 
                   product(xrange(1,6), repeat=3) ): 
    assert all(u == 1)

  iterator = NDimIterator(5, 5, 5)
  u = iterator.next()
  try: u[1] = 2
  except RuntimeError: pass
  else: raise Exception()

  try: NDimIterator(5, 0, 1)
  except ValueError: pass
  else: raise Exception()

  try: NDimIterator(5, 'a', 1)
  except TypeError: pass
  else: raise Exception()

if __name__ == '__main__':
  test_isinteger()
  test_ndimiterator()
