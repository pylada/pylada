//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>

#include "create_clusters.hpp"
#include "../create_clusters.h"


namespace LaDa
{
  namespace Python
  {
    boost::shared_ptr< CE::t_ClusterClasses >
      create_clusters( Crystal::Lattice const &_lat, size_t _max, size_t _order, size_t _site )
        { return CE::create_clusters(_lat, _max, _order, _site ); }
    boost::shared_ptr< CE::t_Clusters >
      create_clusters2( Crystal::Lattice const &_lat, CE::Cluster const &_c, size_t _site )
        { return CE::create_clusters(_lat, _c, _site ); }

    void expose_create_clusters()
    {
      namespace bp = boost::python;

      bp::def
      (
        "create_clusters",
        &create_clusters,
        (
          bp::arg("lattice"),
          bp::arg("nth_shell"),
          bp::arg("order"),
          bp::arg("site") = 0
        ),
        "Returns an array of arrays of equivalent clusters of given order, origin,"
        " and maximum distance (in number of shells) from the origin."
      );

      bp::def
      (
        "create_clusters",
        &create_clusters2,
        (
          bp::arg("lattice"),
          bp::arg("cluster"),
          bp::arg("site") = 0
        ),
        "Returns an array clusters equivalent to the input cluster."
      );
    }

  }
} // namespace LaDa
