def functional():
  from pylada.vff.functional import Functional
  vff = Functional()
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff

def test_inas(epsilon = 1e-4, cartesian=False):
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import exists
  from numpy import dot, sqrt, all, abs, identity
  from pylada.crystal.binary import zinc_blende
  from quantities import angstrom
  from pylada.misc import Changedir

  vff = functional()
  vff.cartesian=cartesian

  directory = '/tmp/test' #mkdtemp()
  if directory == '/tmp/test':
    if exists(directory): rmtree(directory)
    with Changedir(directory) as cwd: pass
  try: 
    structure = zinc_blende()
    structure[0].type = 'In'
    structure[1].type = 'As'
    structure.scale = 2.62332 * 4 / sqrt(3)  * angstrom
    out = vff(structure, outdir=directory, tol=1e-9, overwrite=True)
    startdate = out.start_date
    enddate   = out.end_date
    assert out.success
    assert all(abs(out.structure.cell - structure.cell) < 1e-8)
    assert all(abs(out.structure[0].pos - structure[0].pos) < 1e-8)
    assert all(abs(out.structure[1].pos - structure[1].pos) < 1e-8)
    out = vff(structure, outdir=directory, tol=1e-9, overwrite=False)
    assert out.start_date == startdate
    assert out.end_date == enddate
    assert startdate != enddate
    out = vff(structure, outdir=directory, tol=1e-9, overwrite=True)
    assert out.start_date != startdate

    diff = structure.copy()
    epsilon = identity(3, dtype='float64')                                     \
              + [[0.02, 0.03, -0.1], [0.03, -0.05, -0.05], [-0.1, -0.05, 0.06]]
    epsilon = identity(3, dtype='float64') * 1.1
    diff.cell = dot(epsilon, diff.cell)
    for atom in diff: atom.pos = dot(epsilon, atom.pos)
    newout = vff(diff, outdir=directory, tol=1e-10, overwrite=True)
    assert all(abs(newout.structure.cell - structure.cell) < 1e-8)
    assert all(abs(newout.structure[0].pos - structure[0].pos) < 1e-8)
    assert all(abs(newout.structure[1].pos - structure[1].pos) < 1e-8)
    assert abs(newout.structure.energy) < 1e-8

    diff = structure.copy()
    epsilon = identity(3, dtype='float64')                                     \
              + [[0.02, 0, 0], [0, -0.05, 0], [0, 0, 0.06]]
    diff.cell = dot(epsilon, diff.cell)
    for atom in diff: atom.pos = dot(epsilon, atom.pos)
    newout = vff(diff, outdir=directory, tol=1e-10, overwrite=True)
    assert all(abs(newout.structure.cell - structure.cell) < 1e-8)
    assert all(abs(newout.structure[0].pos - structure[0].pos) < 1e-8)
    assert all(abs(newout.structure[1].pos - structure[1].pos) < 1e-8)
    assert abs(newout.structure.energy) < 1e-8

    diff = structure.copy()
    epsilon = identity(3, dtype='float64')                                     \
              + [[0.02, 0.05, 0], [0.05, -0.05, 0], [0, 0, 0.06]]
    diff.cell = dot(epsilon, diff.cell)
    for atom in diff: atom.pos = dot(epsilon, atom.pos)
    newout = vff(diff, outdir=directory, tol=1e-10, overwrite=True)
    assert all(abs(newout.structure.cell - structure.cell) < 1e-8)
    assert all(abs(newout.structure[0].pos - structure[0].pos) < 1e-8)
    assert all(abs(newout.structure[1].pos - structure[1].pos) < 1e-8)
    assert abs(newout.structure.energy) < 1e-8

    if cartesian: return
    diff = structure.copy()
    diff[1].pos += [0.01, -0.005, 0.03]
    newout = vff(diff, outdir=directory, tol=1e-10, overwrite=True)
    assert all(abs(newout.structure.cell - structure.cell) < 1e-8)
    assert all(abs(newout.structure[1].pos - structure[1].pos) < 1e-8)
    assert all(abs(newout.structure[0].pos - structure[0].pos) < 1e-8)
    assert abs(newout.structure.energy) < 1e-8

    diff = structure.copy()
    epsilon = identity(3, dtype='float64')                                     \
              + [[0.02, 0.05, 0], [0.05, -0.05, 0], [0, 0, 0.06]]
    diff.cell = dot(epsilon, diff.cell)
    for atom in diff: atom.pos = dot(epsilon, atom.pos)
    diff[1].pos += [0.008, -0.005, 0.03]
    newout = vff(diff, outdir=directory, tol=1e-10, overwrite=True)
    assert all(abs(newout.structure.cell - structure.cell) < 1e-8)
    assert all(abs(newout.structure[1].pos - structure[1].pos) < 1e-8)
    assert all(abs(newout.structure[0].pos - structure[0].pos) < 1e-8)
    assert abs(newout.structure.energy) < 1e-8
  finally:
    if directory != '/tmp/test': rmtree(directory)
    
def crystal():
  from pylada.crystal import Structure

  structure = Structure( 10.0, 0.5, 0.5, 
                         0.00, 0.0, 0.5,
                         0.00, 0.5, 0.0, scale=6.5 )\
                       .add_atom(pos=(0.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(0.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(1.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(1.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(2.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(2.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(3.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(3.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(4.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(4.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(5.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(5.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(6.00, 0.00, 0.00), type="In") \
                       .add_atom(pos=(6.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(7.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(7.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(8.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(8.25, 0.25, 0.25), type="As") \
                       .add_atom(pos=(9.00, 0.00, 0.00), type="Ga") \
                       .add_atom(pos=(9.25, 0.25, 0.25), type="As")  
  return structure
def test_ingaas():
  from tempfile import mkdtemp
  from shutil import rmtree
  from os.path import exists
  from pylada.misc import Changedir

  vff = functional()

  directory = '/tmp/test' #mkdtemp()
  if directory == '/tmp/test':
    if exists(directory): rmtree(directory)
    with Changedir(directory) as cwd: pass
  try: 
    structure = crystal()
    vff.method = 'Powell'
    newout = vff(structure, outdir=directory, tol=1e-5, maxiter=200, overwrite=True, verbose=True)
    assert newout.success
    print newout.structure
    print newout.optimize
    print newout.success
  finally:
    if directory != '/tmp/test': rmtree(directory)
    
      
if __name__ == '__main__':
  import scipy
  if tuple(int(u) for u in scipy.__version__.split('.')) >= (0, 11, 0):
    test_inas(1e-7, False)
    test_inas(1e-7, True)
    test_ingaas()
