//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include <opt/debug.h>
#include <opt/tinyxml.h>

#include "electric_dipole.h"


namespace LaDa
{
  namespace GA
  {
    namespace Keepers
    {
      bool OscStrength :: Load( const TiXmlElement &_node )
      {
        if( not _node.Attribute( "oscillator" ) ) 
        {
          std::cerr << "Could not find oscillator attribute.\n";
          return false;
        }
        try
        {
          osc_strength = boost::lexical_cast<types::t_real>
                                            ( _node.Attribute( "oscillator" ) );
        }
        catch( std::exception &_e )
        {
          std::cerr << __SPOT_ERROR << "Error while parsing oscillator strength.\n" << _e.what();
          return false;
        }
        return true;
      }
      bool OscStrength :: Save( TiXmlElement &_node ) const
      {
        _node.SetAttribute( "oscillator", boost::lexical_cast< std::string >( osc_strength ) );
        return true;
      }
    }
    namespace OscStrength
    {
      bool Darwin :: Load( const TiXmlElement &_node )
      {
        degeneracy = types::tolerance;
        const TiXmlElement *parent = opt::find_functional_node( _node, "oscillator strength" );
        if( not parent ) return false;
        if( not parent->Attribute( "degeneracy" ) ) return false;
        degeneracy = boost::lexical_cast< types::t_real >
                                        ( parent->Attribute("degeneracy") );
        return true;
      }
    }
  }
}