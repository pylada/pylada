def functional():
  from pylada.vff.vff import Vff
  vff = Vff()
  vff["In", "As"] = 2.62332, 21.6739, -112.0, 150.0
  vff["Ga", "As"] = 2.44795, 32.1530, -105.0, 150.0
  vff["As", "Ga", "As"] = "tet", -4.099, 9.3703
  vff["Ga", "As", "Ga"] = "tet", -4.099, 9.3703
  vff["In", "As", "In"] = "tet", -5.753, 5.7599
  vff["As", "In", "As"] = "tet", -5.753, 5.7599
  vff["Ga", "As", "In"] = -0.35016, -4.926, 7.5651

  return vff

def test_gradients(epsilon = 1e-4):
  from numpy import abs, dot, array, sqrt
  from pylada.crystal.binary import zinc_blende

  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 6.5 #2.62332 * 2 / sqrt(3)  / 0.529177

  out = vff(structure)
  for atom, check in zip(structure, out):
    oldpos = atom.pos.copy()
    
    for dir in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
      dir = array(dir) / sqrt(dot(dir, dir))
      atom.pos = epsilon * dir + oldpos
      xplus = vff(structure).energy
      atom.pos = -epsilon * dir + oldpos
      xminus = vff(structure).energy
      assert abs(xplus - xminus) < 1e-8


  strain = array([[1e0, 0.1, 0], [0.1, 1e0, 0], [0, 0, 1e0]])
  structure.cell = dot(strain, structure.cell)
  for atom in structure: atom.pos = dot(strain, atom.pos)
  out = vff(structure)
  for atom, check in zip(structure, out):
    oldpos = atom.pos.copy()
    
    for dir in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
      dir = array(dir) / sqrt(dot(dir, dir))
      atom.pos = epsilon * dir + oldpos
      xplus = vff(structure).energy
      atom.pos = -epsilon * dir + oldpos
      xminus = vff(structure).energy

      deriv = (xplus - xminus).magnitude / (2e0*structure.scale * epsilon)
      assert abs(deriv.magnitude - dot(check.gradient, dir)) < 1e2*epsilon

def test_stress(epsilon = 1e-4):
  from numpy import abs, dot, array, identity
  from numpy.linalg import det
  from quantities import angstrom
  from pylada.crystal.binary import zinc_blende

  vff = functional()

  structure = zinc_blende()
  structure[0].type = 'In'
  structure[1].type = 'As'
  structure.scale = 6.5

  out = vff(structure)
  for i in xrange(3):
    for j in xrange(i+1, 3):
      assert abs(out.stress[i, j]) < 1e-8
      assert abs(out.stress[j, i]) < 1e-8

    strain = identity(3, dtype='float64')
    strain[i,i] += epsilon
    newstruct = structure.copy()
    newstruct.cell = dot(strain, newstruct.cell)
    for atom in newstruct: atom.pos = dot(strain, atom.pos)
    xplus = vff(newstruct).energy

    strain = identity(3, dtype='float64')
    strain[i,i] -= epsilon
    newstruct = structure.copy()
    newstruct.cell = dot(strain, newstruct.cell)
    for atom in newstruct: atom.pos = dot(strain, atom.pos)
    xminus = vff(newstruct).energy

    stress = (xplus - xminus) / epsilon * 0.5                                  \
             * -1e0 / det(structure.cell*structure.scale)*angstrom**(-3)
    iival = out.stress[i, i] 
    if not hasattr(iival, 'magnitude'): iival *= out.stress.units 
    assert abs(stress - iival) < 1e2 * epsilon


  ostrain = array([[1e0, 0.1, 0.2], [0.1, 1e0, -0.05], [0.2, -0.05, 1e0]])
  structure.cell = dot(ostrain, structure.cell)
  for atom in structure: atom.pos = dot(ostrain, atom.pos)
  out = vff(structure)
  for i in xrange(3):
    for j in xrange(i, 3):
      strain = identity(3, dtype='float64')
      if i == j: strain[i, i] += epsilon
      else:
        strain[i,j] += 0.5*epsilon
        strain[j,i] += 0.5*epsilon

      newstruct = structure.copy()
      newstruct.cell = dot(strain, newstruct.cell)
      for atom in newstruct: atom.pos = dot(strain, atom.pos)
      xplus = vff(newstruct).energy

      strain = identity(3, dtype='float64')
      if i == j: strain[i, i] -= epsilon
      else:
        strain[i,j] -= 0.5*epsilon
        strain[j,i] -= 0.5*epsilon
      newstruct = structure.copy()
      newstruct.cell = dot(strain, newstruct.cell)
      for atom in newstruct: atom.pos = dot(strain, atom.pos)
      xminus = vff(newstruct).energy

      stress = (xplus - xminus) / epsilon * 0.5                                \
               * -1e0 / det(structure.cell*structure.scale)*angstrom**(-3)
      assert abs(out.stress[i,j]-out.stress[j,i]) < 1e-8
      ijval = out.stress[i, j] 
      if not hasattr(ijval, 'magnitude'): ijval *= out.stress.units 
      assert abs(stress - ijval) < 1e2 * epsilon

      
if __name__ == '__main__':
  test_gradients(1e-9)
  test_stress(1e-8)
