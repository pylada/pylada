#include "LaDaConfig.h"


#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/handle.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include <python/misc.hpp>
#include <python/xml.hpp>

#include "../physics.h"

namespace LaDa
{
  namespace Python
  {
    void expose_physics()
    {
      namespace bp = boost::python;
      bp::def("Z", Physics::Atomic::Z,
              "Given an atomic symbol, returns the atomic number.");
      bp::def("Symbol", Physics::Atomic::Symbol,
              "Given an atomic number, returns the atomic symbol.");
      types::t_unsigned (*ptr_charge) (const std::string &) = &Physics::Atomic::Charge;
      bp::def("Charge", ptr_charge,
              "Given an atomic symbol, returns the number of valence electrons." );
      bp::def("Mass", Physics::Atomic::Mass,
              "Given an atomic symbol, returns the atomic mass.");
      bp::def("a0", Physics::a0,
              "Returns the Bhor radius in A, nm, m, or cm" );
      bp::def("Hartree", Physics::Hartree,
              "Returns the Hartree energy in eV, Rydberg, or Hartree" );
      bp::def("Rydberg", Physics::Rydberg,
              "Returns the Rydberg energy in eV, Rydberg, or Hartree" );
      bp::def("emass", Physics::emass, "Returns the mass of the electron in eV, amu, kg, g, MeV." );
      bp::def("hbar", Physics::hbar, "Returns hbar in  eV*s, erg*s, J*s." );
      bp::def("planck", Physics::Planck, "Returns the Planck constant in  eV*s, erg*s, J*s, Ry, H." );
      bp::def("vacuum_permittivity", Physics::vacuum_permittivity,\
              "Returns the permittivity of vacuum in mks, SI, and au.");
      bp::def("e", Physics::e,\
              "Returns the elementary charge in mks, SI, cgs, and au.");
    }

  }
} // namespace LaDa

BOOST_PYTHON_MODULE(_physics)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "Physics quantities. ";
  scope.attr("__docformat__") = "restructuredtext en";
  bp::docstring_options doc_options(true, false);
  LaDa::Python::expose_physics();
}
