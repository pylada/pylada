#ifndef LADA_CRYSTAL_ATOM_BASE_H
#define LADA_CRYSTAL_ATOM_BASE_H

#include "LaDaConfig.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <complex>

#include <boost/serialization/access.hpp>

#include <opt/types.h>
#include <math/eigen.h>
#include <math/fuzzy.h>

#include <math/serialize.h>
# ifdef LADA_WITH_LNS
#  include <load_n_save/lns.h>
# endif

namespace LaDa
{
  namespace Crystal
  {

    //! \brief Describes an atom.
    //! \details An atom consists of a position and a type. The position should
    //!          always be in cartesian units. The type can be anything, from a
    //!          string with the symbol of the atom, to an double wich codes
    //!          for the atomic type somehow, to a vector of strings which
    //!          describe the possible occupations of the atomic position. To
    //!          this end, the type is a template type \a T_TYPE. 
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation or site ).
    template<class T_TYPE>
    class AtomBase
    {
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
#       friend class load_n_save::access;
#     endif
      public:
        //! The type of the occupation
        typedef T_TYPE t_Type;

      public:
        //! The atomic position in cartesian coordinate.
        math::rVector3d pos;
        //! The atomic occupation.
        t_Type  type;
        
        //! Constructor
        AtomBase() : pos(math::rVector3d(0,0,0)), type() {};
        //! Constructor and Initializer
        explicit AtomBase   ( const math::rVector3d &_pos, t_Type _type) 
                          : pos(_pos), type(_type) {};
        //! Copy Constructor
        AtomBase(const AtomBase &_c) : pos(_c.pos), type(_c.type) {};

      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { _ar & pos; _ar & type; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#       endif
    };

#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_TYPE> template<class T_ARCHIVE>
        bool AtomBase<T_TYPE> lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          namespace lns = LaDa :: load_n_save;
          return _ar & 
                 (
                   lns::section("Atom")
                     << lns::option( "pos", lns::tag=lns::required, lns::action=pos,
                                     lns::help="Cartesian position in Anstrom." )
//                    << lns::option( "type", lns::tag=lns::required, lns::action=type,
//                                    lns::help="Atomic specie, or any string." )
                 );
         }
#   endif

  } // namespace Crystal
} // namespace LaDa
  
#endif