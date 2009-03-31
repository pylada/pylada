//
//  Version: $Id$
//
#ifndef _LADA_ATAT_SERIALIZE_H_
#define _LADA_ATAT_SERIALIZE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vectmac.h"

namespace boost {
  namespace serialization {

    //! Serializes atat real vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::rVector3d & g, const unsigned int version)
     { ar & g.x; }
     //! Serializes atat integer vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::iVector3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat real matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::rMatrix3d & g, const unsigned int version)
     { ar & g.x; }
    //! Serializes atat integer matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::atat::iMatrix3d & g, const unsigned int version)
     { ar & g.x; }
  }

}

#endif