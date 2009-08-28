//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>

#include "values.h"
#include "../collapse/values.h"
#include "../sum_of_separables.h"
#include "../representation.h"
#include "../collapse/fitting_set.h"
#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>


namespace LaDa
{
  namespace Python
  {
    void expose_values()
    {
      namespace bp = boost::python;
      typedef LaDa::atomic_potential::collapse::Values Values;
      
      expose_vector< std::vector< std::vector<double> > >
         ("Values_vec_vec", "Implementation only.");
      expose_vector<Values::t_CoordRankValues::value_type>
         ("Values_vec_vec_vec", "Implementation only.");
      expose_vector<Values::t_CoordRankValues>
        ("Values_vec_vec_vec_vec", "Implementation only.");
      expose_vector<Values::t_FunctionValues>
        ("Values_vec_vec_vec_vec_vec", "Implementation only.");
 
      bp::class_<Values>("Values", "Holds fitting related info. For Debugging only.")
        .add_property("coord_rank_values", &Values::coord_rank_values_)
        .add_property("rank_values", &Values::rank_values_)
        .add_property("function_values", &Values::function_values_)
        .def("add", &Values::add);
    }

  }
} // namespace LaDa
