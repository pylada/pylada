def functional():
  from pylada.vff.functional import Functional
  vff = Functional()
  vff["Si", "Si"] = 2.35905923, 45.375785351, -150.363195485, 519.230157133
  vff["Ge", "Ge"] = 2.44167023, 37.816422666, -124.830189020, 250.179861133
  vff["Si", "Ge"] = 2.39642606, 41.875535703, -145.994430120, 366.853558523
  vff["Si", "Si", "Si"] = "tet", 0, 12.579328566, 44.930324684, 359.422663897
  vff["Ge", "Ge", "Ge"] = "tet", 0, 10.383093493, 55.677487481
  vff["Si", "Ge", "Si"] = "tet", 0, 11.902809727, 20.193404352
  vff["Ge", "Si", "Ge"] = "tet", 0, 11.902809727, 20.193404352
  vff["Si", "Si", "Ge"] = "tet", 0, 12.245971457, 32.561864518, 179.706525419
  vff["Si", "Ge", "Ge"] = "tet", 0, 11.147853920, 37.930591733

  return vff

def crystal():
  from pylada.crystal import Structure

  structure = Structure( 10.0, 0.5, 0.5, 
                         0.00, 0.0, 0.5,
                         0.00, 0.5, 0.0, scale=5.5 )\
                       .add_atom(pos=(0.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(0.25, 0.25, 0.25), type="Ge") \
                       .add_atom(pos=(1.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(1.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(2.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(2.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(3.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(3.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(4.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(4.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(5.00, 0.00, 0.00), type="Si") \
                       .add_atom(pos=(5.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(6.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(6.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(7.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(7.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(8.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(8.25, 0.25, 0.25), type="Si") \
                       .add_atom(pos=(9.00, 0.00, 0.00), type="Ge") \
                       .add_atom(pos=(9.25, 0.25, 0.25), type="Si")  
  return structure


def test_relaxepi(direction, cartesian=True, scale=5.45):
  from tempfile import mkdtemp
  from os.path import exists
  from shutil import rmtree
  from numpy import array, outer, abs, all, dot, cross, max
  from numpy.linalg import norm
  from pylada.misc import Changedir

  direction = array(direction) / norm(direction)
  vff = functional()
  vff.cartesian = cartesian
  vff.direction = direction
  structure = crystal()
  structure.scale = scale
  template = outer(direction, direction)

  def rmdir(cell): return cell - dot(template.T, dot(cell, template))
  
  directory = '/tmp/test/' # mkdtmp()
  if directory == '/tmp/test/':
    if exists(directory): rmtree(directory)
    with Changedir(directory) as cwd: pass
  try: 
    out = vff(structure, outdir=directory, overwrite=True, verbose=True)
    assert out.success
    diff = out.structure.cell - out.input_structure.cell
    maxi = max(abs(diff))
    assert all(abs(cross(direction, diff[:, 0])) < maxi * 1e-3)
    assert all(abs(cross(direction, diff[:, 1])) < maxi * 1e-3)
    assert all(abs(cross(direction, diff[:, 2])) < maxi * 1e-3)
  finally:
    if directory != '/tmp/test/': rmtree(directory)

if __name__ == '__main__':
  import scipy
  if tuple(int(u) for u in scipy.__version__.split('.')) >= (0, 11, 0):
    test_relaxepi([0, 0, 1], True, 5.35)
    test_relaxepi([0, 1, 1], True, 5.35)
    test_relaxepi([1, 1, 1], True, 5.35)
