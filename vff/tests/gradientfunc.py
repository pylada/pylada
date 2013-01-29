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

def test_relaxall(cartesian):
  """ Tests gradients as implemented in the functional. """
  from numpy import abs
  from pylada.vff.functional import Functional

  vff = functional()
  vff.cartesian = cartesian
  structure = crystal()
  epsilon = 1e-5
  
  # gets functions for current structure
  x, jacobian, energy, make_structure                                          \
       = Functional(copy=vff)._getfuncs_relaxall(structure)
  stress, forces = vff.jacobian(structure)

  jac0 = jacobian(x)
  for i in xrange(6):
    xminus = x.copy()
    xminus[i] = xminus[i] - epsilon * 0.5
    eminus = energy(xminus)
    xplus = x.copy()
    xplus[i] = xplus[i] + epsilon * 0.5
    eplus = energy(xplus)

    assert abs(jac0[i] - (eplus - eminus) / epsilon ) < 1e-6

  epsilon = 1e-7
  for i in xrange(6, len(x)):
    xminus = x.copy()
    xminus[i] = xminus[i] - epsilon * 0.5
    eminus = energy(xminus)
    xplus = x.copy()
    xplus[i] = xplus[i] + epsilon * 0.5
    eplus = energy(xplus)

    assert abs(jac0[i] - (eplus - eminus) / epsilon ) < 1e-6

if __name__ == '__main__':
  test_relaxall(False)
  test_relaxall(True)