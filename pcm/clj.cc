#include "LaDaConfig.h"

#include <fstream>
#include <string>
#include <algorithm>

#include <boost/spirit/include/classic_primitives.hpp>
#include <boost/spirit/include/classic_numerics.hpp>
#include <boost/spirit/include/classic_actions.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_operators.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lambda/core.hpp>

#include <physics/physics.h>

#include "clj.h"

namespace LaDa
{
  namespace Models
  {
    bool Clj :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Clj" );
      if( not parent ) return false;
      if( not Ewald::Load( *parent ) ) return false;
      if( not LennardJones :: Load( *parent ) ) return false;
      return true;
    }

    Clj :: t_Return Clj :: call_(const t_Arg& _in, t_Arg& _out, size_t _which) const
    {
      _out.cell = math::rMatrix3d::Zero();
      foreach( t_Arg :: t_Atom &atom, _out.atoms )
        atom.pos = math::rVector3d(0,0,0);

      _out.energy = 0e0;
      if( _which & 1 ) _out.energy += LennardJones::operator()( _in, _out );
      if( _which & 2 ) _out.energy += Ewald::operator()( _in, _out );
      return _out.energy;
    }

    Clj :: t_Return Clj :: gradient(const t_Arg& _in, t_Arg& _out) const
    {
      operator()( _in, _out );
      _out.cell = -_out.cell * (~(!_in.cell));
      t_Arg :: t_Atoms :: iterator i_force = _out.atoms.begin();
      const t_Arg :: t_Atoms :: iterator i_force_end = _out.atoms.end();
      const math::rMatrix3d matrix( -_in.scale * (~_in.cell) );
      for(; i_force != i_force_end; ++i_force )
        i_force->pos = matrix * i_force->pos;
      return _out.energy;
    }

    namespace details
    {
      // contains atom specific quantities.
      struct atomic_specie
      {
        size_t atomic_number;
        types::t_real radius;
        types::t_real charge;
      };
      bool operator==( const atomic_specie &_a, const atomic_specie &_b )
        { return _a.atomic_number == _b.atomic_number; } 
      bool operator==( const atomic_specie &_a, const size_t &_b )
        { return _a.atomic_number == _b; } 
     
      // contains bond specific quantities.
      struct bond_type
      {
        size_t a;
        size_t b;
        types::t_real epsilon;
        types::t_real rsigma;
      };
      bool operator==( const bond_type &_a, const bond_type &_b )
      {
        return (_a.a == _b.a and  _a.b == _b.b ) or (_a.a == _b.b and  _a.b == _b.a );
      }
    }

    void read_fortran_input( Clj &_clj, std::vector<std::string> &_atoms, 
                             const boost::filesystem::path &_path )
    {
      __TRYBEGIN

      namespace bsc = boost::spirit::classic;
      namespace fs = boost::filesystem;  
      __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                  _path << " is neither a regulare file nor a system link.\n" )
      std::ifstream file( _path.string().c_str(), std::ifstream::in );
      std::string line;

      __DOASSERT( file.bad(), "Could not open file " + _path.string() + "\n" );
      
      // Number of atomic species.
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      size_t nbspecies;
      __DOASSERT
      (
        not bsc::parse
        (
          line.c_str(),
           bsc::uint_p[ bsc::assign_a( nbspecies ) ],
           bsc::space_p
        ).hit,
        "Could not parse number of species.\n" 
      )

      // Atomic types
      __DOASSERT( nbspecies == 0,
                  "Number of species set to zero in input " + _path.string() + ".\n" )
      std::vector< details::atomic_specie > species;
      _atoms.clear();
      for( size_t i(0); i < nbspecies; ++i )
      {
        details::atomic_specie atom;
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
          
                bsc::uint_p[ bsc::assign_a( atom.atomic_number ) ] 
             >> bsc::ureal_p[ bsc::assign_a( atom.radius ) ] 
             >> bsc::real_p[ bsc::assign_a( atom.charge ) ],
             bsc::space_p
          ).hit,
          "Could not parse atomic type.\n" 
        )
        __DOASSERT( species.end() != std::find( species.begin(), species.end(), atom ),
                    "Same specie specified twice.\n" )
        species.push_back( atom );
        __DOASSERT( Physics::Atomic::Symbol(atom.atomic_number) == "error",
                    "Unknown atom with z=" << atom.atomic_number << ". Please complain.\n" )
        _atoms.push_back( Physics::Atomic::Symbol( atom.atomic_number ) );
      }
      
      // Bond types
      __DOASSERT( nbspecies == 0,
                  "Number of species set to zero in input " + _path.string() + ".\n" )
      std::vector< details::bond_type > bonds;
      for( size_t i(0); i < nbspecies * nbspecies; ++i )
      {
        details::bond_type bond;
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
          
                bsc::uint_p[ bsc::assign_a( bond.a ) ] 
             >> bsc::uint_p[ bsc::assign_a( bond.b ) ] 
             >> bsc::ureal_p[ bsc::assign_a( bond.epsilon ) ] 
             >> bsc::real_p[ bsc::assign_a( bond.rsigma ) ],
             bsc::space_p
          ).hit,
          "Could not parse atomic type.\n" 
        )
        const std::vector< details::bond_type > :: const_iterator i_found
        ( 
          std::find( bonds.begin(), bonds.end(), bond )
        );
        if( i_found != bonds.end() )
        {
          __DOASSERT( i_found->a == bond.a and i_found->b == bond.b,
                      "Same bond specified twice (bond " << i << ").\n") 
          __DOASSERT( std::abs( i_found->epsilon - bond.epsilon ) > 1e-12,
                      "Same bond, different epsilon (bond " << i << ").\n") 
          __DOASSERT( std::abs( i_found->rsigma - bond.rsigma ) > 1e-12,
                      "Same bond, different rsigma (bond " << i << ").\n") 
          continue;
        }
        bonds.push_back( bond );
      }

      // lj real space cutoff
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      __DOASSERT
      (
        not bsc::parse
        (
          line.c_str(),
           bsc::ureal_p[ bsc::assign_a( _clj.LennardJones::rcut_ ) ],
           bsc::space_p
        ).hit,
        "Could not parse real-space cut-off for lj.\n" 
      )
      // lj real space cutoff
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      __DOASSERT
      (
        not bsc::parse
        (
          line.c_str(),
           bsc::ureal_p[ bsc::assign_a( _clj.Ewald::cutoff_ ) ],
           bsc::space_p
        ).hit,
        "Could not parse real-space cut-off for ewald.\n" 
      )

      // PEGS.
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      types::t_real pegs;
      __DOASSERT
      (
        not bsc::parse
        (
          line.c_str(),
           bsc::ureal_p[ bsc::assign_a( pegs ) ],
           bsc::space_p
        ).hit,
        "Could not parse real-space cut-off for ewald.\n" 
      )

      // Now converts to Ewald.
      _clj.Ewald::charges.clear(); 
      foreach( const details::atomic_specie & specie, species )
      {
        __DOASSERT( Physics::Atomic::Symbol( specie.atomic_number ) == "error",
                    "Unknow atomic number " << specie.atomic_number << ", please complain.\n" )
        _clj.Ewald::charges[ Physics::Atomic::Symbol( specie.atomic_number ) ] = specie.charge;
      }

      // Now converts to LJ.
      _clj.LennardJones::bonds.clear();
      foreach( const details::bond_type & bond, bonds )
      {
        __DOASSERT( Physics::Atomic::Symbol( bond.a ) == "error",
                    "Unknow atomic number " << bond.a << ", please complain.\n" )
        __DOASSERT( Physics::Atomic::Symbol( bond.b ) == "error",
                    "Unknow atomic number " << bond.b << ", please complain.\n" )
        __DOASSERT( species.end() == find( species.begin(), species.end(), bond.a ), 
                    "Atomic type " << bond.a << " incomplete in input.\n" )
        __DOASSERT( species.end() == find( species.begin(), species.end(), bond.b ), 
                    "Atomic type " << bond.a << " incomplete in input.\n" )

        const std::string bondtype
        ( 
          _clj.LennardJones::bondname
          ( 
             Physics::Atomic::Symbol( bond.b ), 
             Physics::Atomic::Symbol( bond.a ) 
          )
        );
        __ASSERT( _clj.LennardJones::bonds.end() != _clj.LennardJones::bonds.find( bondtype ),
                  "Bond already exists.\n" )
        const details::atomic_specie A( *find( species.begin(), species.end(), bond.a ) );
        const details::atomic_specie B( *find( species.begin(), species.end(), bond.b ) );
        const types::t_real radius( ( A.radius + B.radius ) * bond.rsigma );
        const types::t_real radius2( radius * radius );
        const types::t_real radius6( radius2 * radius2 * radius2 );
        const types::t_real radius12( radius6 * radius6 );
        _clj.LennardJones::bonds[ bondtype ].hard_sphere = radius12 * bond.epsilon * 4e0;
        _clj.LennardJones::bonds[ bondtype ].van_der_walls = radius6 * bond.epsilon * (1e0 - pegs) * 4e0;
      }

      // Hard coded in fortran.
      boost::tuples::get<0>(_clj.LennardJones::mesh_) = 3;
      boost::tuples::get<1>(_clj.LennardJones::mesh_) = 3;
      boost::tuples::get<2>(_clj.LennardJones::mesh_) = 3;

      __TRYEND(, "Could not parse input.\n" )
    }

    void Clj :: check_coherency() const
    {
      namespace bl = boost::lambda;
      LennardJones :: check_coherency();
      std::vector< std::string > ewald;
      Ewald :: t_Charges :: const_iterator i_ewald = Ewald :: charges.begin();
      Ewald :: t_Charges :: const_iterator i_ewald_end = Ewald :: charges.end();
      for(; i_ewald != i_ewald_end; ++i_ewald) 
        ewald.push_back( i_ewald->first );

      std::vector< std::string > lj;
      LennardJones :: t_Bonds :: const_iterator i_bond = LennardJones :: bonds.begin();
      LennardJones :: t_Bonds :: const_iterator i_bond_end = LennardJones :: bonds.end();
      for(; i_bond != i_bond_end; ++i_bond) 
      {
        const std::string A( LennardJones :: extract_atomA( i_bond->first ) );
        const std::string B( LennardJones :: extract_atomB( i_bond->first ) );
        if( lj.end() == std::find( lj.begin(), lj.end(), A ) ) lj.push_back( A );
        if( lj.end() == std::find( lj.begin(), lj.end(), B ) ) lj.push_back( B );
      }

      bool result = true;
      std::sort( ewald.begin(), ewald.end() );
      std::sort( lj.begin(), lj.end() );
      std::vector<std::string> diff;
      std::set_difference( ewald.begin(), ewald.end(), lj.begin(), lj.end(),
                           std::back_inserter( diff ) );
      if( diff.size() != 0 )
      {
        result = false;
        std::cerr << "These atomic types are present in Ewald sum, but not in LennardJones.\n";
        std::for_each( diff.begin(), diff.end(), std::cerr << bl::_1 << bl::constant( " " ) );
        std::cerr << "\n";
      }
      diff.clear();
      std::set_difference( lj.begin(), lj.end(), ewald.begin(), ewald.end(), 
                           std::back_inserter( diff ) );
      if( diff.size() != 0 )
      {
        result = false;
        std::cerr << "These atomic types are present in LennardJones sum, "
                     "but not in Ewald sum.\n";
        std::for_each( diff.begin(), diff.end(), std::cerr << bl::_1 << bl::constant( " " ) );
        std::cerr << "\n";
      }
      __DOASSERT( not result, "" );
    }

    std::ostream& operator<<( std::ostream& _stream, const Clj& _clj )
    {
      _stream << "Lennard-Jones + Coulomb functional:\n  LennardJones Parameters:\n";
      Clj::LennardJones::t_Bonds :: const_iterator i_bond = _clj.bonds.begin();
      Clj::LennardJones::t_Bonds :: const_iterator i_bond_end = _clj.bonds.end();
      for(; i_bond != i_bond_end; ++i_bond )
        _stream << "    "
                << i_bond->first << ": " << i_bond->second.hard_sphere << "/r^12 - " 
                << i_bond->second.van_der_walls << "/r^6\n";
      Clj::Ewald::t_Charges :: const_iterator i_charge = _clj.charges.begin();
      Clj::Ewald::t_Charges :: const_iterator i_charge_end = _clj.charges.end();
      _stream << "  Coulomb terms:\n";
      for(; i_charge != i_charge_end; ++i_charge )
        _stream << "    "
                << i_charge->first << ": " << i_charge->second << "/r^2\n";
      return _stream;
    }
  } // namespace CLJ
} // namespace LaDa