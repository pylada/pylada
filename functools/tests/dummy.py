class Base(object):
  def __init__(self, a=2, b='b'):
    super(Base, self).__init__()
    self.a = a
    self.b = b
  def __call__(self, structure, outdir=None):
    self.b += 1

iterator = -1

def iter_func(self, structure, outdir=None, other=True):
  assert other == False
  global iterator
  for iterator in xrange(1, 5): 
    self.a *= iterator
    self(structure, outdir)
    yield self
def func(self, structure, outdir=None, other=True):
  assert other == False
  self(structure, outdir)
  a = self.a
  for i, f in enumerate(iter_func(self, structure, outdir, other=other)): 
    assert f.a / a == i + 1
    a = f.a
  self.a = 0
  return not structure
    


