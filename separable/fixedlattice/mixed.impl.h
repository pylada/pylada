//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/lambda/bind.hpp>

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE) || defined(INCOLLAPSE2)
#   error "Macros with same names."
# endif
# define COLHEAD \
    MixedApproach<T_TRAITS> 
#  define INCOLLAPSE( var ) \
     template< class T_TRAITS >  var COLHEAD
#  define INCOLLAPSE2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 COLHEAD

    INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
      :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        __ASSERT( not t_ColBase::separables_, "Function pointer not set.\n" )
        // Loop over inequivalent configurations.
        t_Vector X( dof() );
        std::fill( _A.data().begin(), _A.data().end(), 
                   typename t_Matrix::value_type(0) );
        std::fill( _b.data().begin(), _b.data().end(), 
                   typename t_Vector::value_type(0) );
        const bblas::range colrange( 0, t_ColBase::dof() );
        typename t_CEBase :: t_Pis :: const_iterator i_pis;
        if( t_CEBase::dof() ) 
        {
          i_pis = t_CEBase::pis.begin();
          __ASSERT( t_CEBase::pis.size() != t_ColBase::mapping().size(), 
                    "Inconsistent sizes.\n" )
        }
        for( size_t i(0); i < t_ColBase::mapping().size(); ++i )
        {
          // allows leave-one-out, or leave-many-out.
          if( t_ColBase::mapping().do_skip(i) ) continue;
      
          // create the X vector.
          bblas::vector_range< t_Vector > colX( X, colrange );
          std::fill( colX.begin(), colX.end(), typename t_Vector::value_type(0) );
          if( t_ColBase::dof() ) t_ColBase::create_X( i, colX );
          if( t_CEBase::dof() )
          {
            std::copy( i_pis->begin(), i_pis->end(), X.begin() + t_ColBase::dof() );
            ++i_pis;
          }
          
      
          _A += t_ColBase::mapping().weight(i) * bblas::outer_prod( X, X ); 
          _b += t_ColBase::mapping().weight(i) * t_ColBase::mapping().target(i) * X;
        }
        __DEBUGTRYEND(, "Error in MixedApproach::create_A_n_b()\n" )
      }

    INCOLLAPSE2( template< class T_MATRIX, class T_VECTOR > void )
      :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        t_ColBase::dim = _dim;
        create_A_n_b( _A, _b );

        const bool dosep( t_ColBase::dof() > 0 );
        if( t_ColBase::dof() )
        {
          const bblas::range seprange( 0, t_ColBase::dof() );
          bblas::matrix_range<T_MATRIX> sepmat( _A, seprange, seprange );
          bblas::vector_range<T_VECTOR> sepvec( _b, seprange );
          t_ColBase::regularization()( _A, _b, _dim); 
        }
        if( t_CEBase::dof() )
        {
          const bblas::range cerange( t_ColBase::dof(), dof() );
          bblas::matrix_range<T_MATRIX> cemat( _A, cerange, cerange );
          bblas::vector_range<T_VECTOR> cevec( _b, cerange );
          t_CEBase :: other_A_n_b( cemat, cevec );

          if( t_ColBase::dof() )
          {
            bblas::matrix_column< t_Matrix > columnd( coefficients(), _dim );
            const bblas::matrix_column< const t_Matrix > column0( coefficients(), 0 );
            std::copy( column0.begin() + t_ColBase::dof(), column0.end(), 
                       columnd.begin() + t_ColBase::dof() );
          }
        }
        __DEBUGTRYEND(, "Error in MixedApproach::operator()()\n" )
      }
    
    INCOLLAPSE( void ) :: randomize( typename t_Vector :: value_type _howrandom )
    {
      namespace bblas = boost::numeric::ublas;
      separables().randomize( _howrandom );
      typedef typename t_Matrix :: value_type t_Type;
      bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
      typename bblas::matrix_column< t_Matrix > :: iterator i_c = column0.begin();
      typename bblas::matrix_column< t_Matrix > :: iterator i_c_end = column0.end();
      for( i_c += t_ColBase::dof(); i_c != i_c_end; ++i_c )
        *i_c = t_Type( opt::random::rng() - 5e-1 ) * _howrandom;
    }


    INCOLLAPSE( void ) :: update_all()
    {
      __DEBUGTRYBEGIN
      if( t_CEBase::dof() )
      {
        namespace bblas = boost::numeric::ublas;
        const bblas::matrix_column< const t_Matrix > column0( coefficients(), 0 );
        for( size_t i(1); i < coefficients().size1(); ++i )
        {
          bblas::matrix_column< t_Matrix > columnd( coefficients(), t_ColBase::dim );
          std::copy( column0.begin() + t_ColBase::dof(), column0.end(), 
                     columnd.begin() + t_ColBase::dof() );
        }
      }
      if( t_ColBase::dof() ) t_ColBase::update_all();
      __DEBUGTRYEND(, "Error in MixedApproach::update_all()\n" )
    }
    INCOLLAPSE( void ) :: update( types::t_unsigned _d )
    {
      __DEBUGTRYBEGIN
      if( t_CEBase::dof() )
      {
        namespace bblas = boost::numeric::ublas;
        const bblas::matrix_column< t_Matrix > columnd( coefficients(),  
                                                        t_ColBase::dim );
        bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
        std::copy( columnd.begin() + t_ColBase::dof(), columnd.end(), 
                   column0.begin() + t_ColBase::dof() );
      }
      if( t_ColBase::dof() ) t_ColBase::update( _d );
      __DEBUGTRYEND(, "Error in MixedApproach::update()\n" )
    }

    INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const 
    {
      __DEBUGTRYBEGIN
      __ASSERT(    coefficients().size1() != dof() 
                or coefficients().size2() != dimensions(),
                "Inconsistent sizes.\n" )
      opt::ErrorTuple error;
      for(size_t n(0); n < t_ColBase::mapping().size(); ++n )
        if( not t_ColBase::mapping().do_skip(n) )
          error += opt::ErrorTuple(   t_ColBase::mapping().target( n )
                                    - evaluate(n),
                                    t_ColBase::mapping().weight( n ) );
      return error;
      __DEBUGTRYEND(, "Error in MixedApproach::evaluate()\n" )
    }

    INCOLLAPSE( typename COLHEAD::t_Matrix::value_type ) :: evaluate(size_t _n) const 
    {
      namespace bblas = boost::numeric::ublas;
      __DEBUGTRYBEGIN
      __ASSERT(    coefficients().size1() != dof() 
                or coefficients().size2() != dimensions(),
                "Inconsistent sizes.\n" )
      typename t_Matrix :: value_type intermed(0);
      // Adds Separable part.
      if( t_ColBase::dof() ) intermed += t_ColBase::evaluate(_n);
      // Adds CE part.
      if( t_CEBase::dof() )
      {
        typedef const bblas::matrix_column<const t_Matrix>  t_Column;
        typedef const bblas::vector_range< t_Column > t_Ecis;
        const bblas::matrix_column< const t_Matrix> col( coefficients(), 0 );
        t_Ecis ecis( col, bblas::range( t_ColBase::dof(), dof() ) );
        intermed += ( bblas::inner_prod( ecis, t_CEBase::pis[ _n ] ) );
      }

      return intermed;
      __DEBUGTRYEND(, "Error in MixedApproach::evaluate()\n" )
    }

    INCOLLAPSE( void ) ::  init( size_t _ranks, size_t _dims )
    {
      __DEBUGTRYBEGIN
      CEFit().init( clusters() );
      if ( not ( _ranks and _dims ) )
      { 
        _ranks = 0;
        _dims = 1;
      }
      namespace bblas = boost::numeric::ublas;
      separables().norms.resize( _ranks );
      coefficients().resize(   clusters_->size() 
                             + _ranks * t_Separables::t_Mapping::D, _dims );
      const bblas::range rows( 0, _ranks * t_Separables::t_Mapping::D );
      const bblas::range columns( 0, _dims );
      separables().coefficients_interface().set( coefficients(), rows, columns );
      t_ColBase :: init( separables_ );
      __DEBUGTRYEND(, "Error in MixedApproach::init()\n" )
    }

    INCOLLAPSE( void ) ::  reassign()
    {
      if( not t_CEBase::dof() ) return;
      namespace bl = boost::lambda;
      __DEBUGTRYBEGIN
      __ASSERT( coefficients().size1() - t_ColBase::dof() != clusters_->size(),
                "Inconsistent sizes.\n")

      typename t_Matrix :: const_iterator1 i_eci = coefficients().begin1() + t_ColBase::dof();
      foreach( t_Clusters::value_type & _clusters, clusters() )
      {
        foreach( ::CE::Cluster & _cluster, _clusters ) _cluster.eci = *i_eci;
        ++i_eci;
      }
      __DEBUGTRYEND(,"Error in MixedApproach::reassign().\n" )
    }

#  undef COLHEAD
#  undef INCOLLAPSE
#  undef INCOLLAPSE2

  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream,
                            const MixedApproach<T_TRAITS> &_col )
  {
    _stream << "Mixed-Approach function: " << _col.separables() << "\n";
    typedef typename MixedApproach<T_TRAITS> :: t_Clusters :: const_iterator t_cit;
 
    t_cit i_clusters = _col.clusters().begin();
    t_cit i_clusters_end = _col.clusters().end();
    for(; i_clusters != i_clusters_end; ++i_clusters )
      _stream << i_clusters->front();
    return _stream;
  }

  template< class T_POSITIONS, class T_CLUSTERS >
    void remove_contained_clusters( const T_POSITIONS &_positions, T_CLUSTERS &_clusters )
    {
      typename T_CLUSTERS :: iterator i_clusters = _clusters.begin();
      typename T_CLUSTERS :: iterator i_clusters_end = _clusters.end();
      for(; i_clusters != i_clusters_end; ++i_clusters )
      {
        if( i_clusters->front().size() < 2 ) continue;
        typename T_CLUSTERS::value_type::iterator i_cluster =  i_clusters->begin();
        typename T_CLUSTERS::value_type::iterator i_cluster_end =  i_clusters->end();
        for(; i_cluster != i_cluster_end; ++i_cluster )
        {
          size_t n(0);
          for(; n < size_t(i_cluster->size()); ++n ) 
          {
            const atat::rVector3d &pos( (*i_cluster)[n] );
            if( Fuzzy::is_zero( atat::norm2( pos ) ) ) continue;
            typename T_POSITIONS :: const_iterator i_pos = _positions.begin();
            typename T_POSITIONS :: const_iterator i_pos_end = _positions.end();
            for(; i_pos != i_pos_end; ++i_pos )
              if( Fuzzy::is_zero( atat::norm2( *i_pos - pos ) ) ) break;
            if( i_pos == i_pos_end ) break;
          }
          if( n == size_t( i_cluster->size() ) ) break;
        }
        if( i_cluster == i_cluster_end ) continue;
        typename T_CLUSTERS :: iterator dummy = i_clusters;
        --i_clusters;
        _clusters.erase( dummy );
        i_clusters_end = _clusters.end();
      }
    }
} // end of CE namespace.
