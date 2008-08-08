//
//  Version: $Id$
//

#include <opt/random.h>
#include <opt/algorithms.h>

namespace Separable
{
# ifdef __DOHALFHALF__
    template< class T_FUNCTION > template< class T_VECTORS >
      void EquivCollapse<T_FUNCTION> :: create_coefs( T_VECTORS& _coefs,
                                                      t_Type _howrandom ) const
      {
        // Initializes norms to 1
        std::fill( function.coefs.begin(), function.coefs.end(), t_Type(1) );
       
        __DEBUGTRYBEGIN
        typedef T_VECTORS t_Vectors;
        typedef typename t_Vectors :: iterator t_solit;
        _coefs.resize( sizes.size() );
        t_solit i_dim = _coefs.begin();
        t_solit i_dim_end = _coefs.end();
        typename t_Base::t_Sizes :: const_iterator i_sizes = sizes.begin();
        types::t_real ndim( sizes.size() );
       
        // loop over dimensions
        for(; i_dim != i_dim_end; ++i_dim, ++i_sizes)
        {
          types::t_int ri = std::accumulate( i_sizes->begin(), i_sizes->end(), 0 );
          i_dim->resize( ri );
          typename t_Base::t_Sizes :: value_type :: const_iterator
            i_size = i_sizes->begin();
          typename t_Base::t_Sizes :: value_type :: const_iterator
            i_size_end = i_sizes->end();
          typename t_Vectors :: value_type :: iterator i_coef = i_dim->begin();
          // loop over ranks.
          for(size_t r(0); i_size != i_size_end; ++i_size, ++r )
          {
            // First creates a normalized random vector.
            namespace bl = boost::lambda;
            typename t_Vectors :: value_type :: iterator i_c( i_coef ); 
            t_Type norm;
            do
            {
              norm = t_Type(0);
              const t_Type range(_howrandom);
              bool centered = true;
              std::for_each 
              (
                i_c, i_c + *i_size, 
                bl::var( norm ) += (
                                     bl::_1 =   bl::constant(range) 
                                              * bl::bind( &opt::random::rng )
                                                + bl::if_then_else_return
                                                  ( 
                                                    bl::var(centered), 
                                                    bl::constant( t_Type(-0.5*range) ),
                                                    bl::constant( t_Type(1e0-0.5*range) )
                                                  ),
                                     bl::var(centered) = !bl::var(centered),
                                     bl::_1 * bl::_1
                                   )
              );
            }
            while( Fuzzy::eq( norm, t_Type(0) ) );
            norm = std::sqrt(norm);
            function[r] = norm;
            norm = t_Type(1) / norm;
            std::for_each 
            (
              i_c, i_c + *i_size, 
              bl::_1 *= bl::constant(norm)
            );
            i_coef += *i_size;
            // Then adds random norm to coefficient
          } // end of loop over ranks.
        } // end of loop over dimensions.
        __DEBUGTRYEND(, "Error while assigning solution coefs to function.\n" )
    }
  #endif

  template<class T_FUNCTION>
    template< class T_VECTORS, class T_VECTOR, class T_EWEIGHTS >
      void EquivCollapse<T_FUNCTION> :: init( const T_VECTORS &_x, const T_VECTOR &_w,
                                              const T_EWEIGHTS &_eweights ) 
      {
        __DEBUGTRYBEGIN
  
        t_Base :: init( _x, _w );
  
        __ASSERT( _w.size() != _eweights.size(),
                  "Inconsistent array-sizes on input.\n" )
        eweights.resize( _eweights.size() );
        typename T_EWEIGHTS :: const_iterator _i_eweights = _eweights.begin();
        t_eWeights :: iterator i_eweights = eweights.begin();
        t_eWeights :: iterator i_eweights_end = eweights.end();
        for(; i_eweights != i_eweights_end; ++i_eweights, ++_i_eweights )
        {
          i_eweights->resize( _i_eweights->size() );
          std::copy( _i_eweights->begin(), _i_eweights->end(), i_eweights->begin() );
        }
  
        __DEBUGTRYEND(, "Error while initializing Collapse "
                        "function for equivalent structures.\n" )
      } // end of init member function.

  template<class T_FUNCTION>  
    template< class T_VECTOR, class T_MATRIX, class T_VECTORS >
      void EquivCollapse<T_FUNCTION>::operator()( T_VECTOR &_b, T_MATRIX &_A, 
                                                  types::t_unsigned _dim,
                                                  const T_VECTOR &_targets,
                                                  const T_VECTORS &_coefs     )
      {
        __DEBUGTRYBEGIN

        __ASSERT( _dim  >= D, "Input dimension out of range.\n" )
        __ASSERT( _coefs.size() != D,
                  "Inconsistent number of dimensions/coefficients.\n" )
        __ASSERT( _targets.size() != weights.size(),
                  "Unexpected array-size on input.\n" )
        __ASSERT( eweights.size() != weights.size(),
                  "Unexpected array-size on input.\n" )
        __ASSERT
        ( 
          _coefs[_dim].size() != std::accumulate( sizes[_dim].begin(),
                                                  sizes[_dim].end(), 0 ),
          "Inconsistent number of ranks/dimensions/coefficients.\n" 
        )
        
        __ASSERT( sizes.size()  != D, "Inconsistent sizes.\n" )
        __ASSERT( factors.size() == _coefs[_dim].size(), "Inconsistent sizes.\n" )
        __ASSERT( expanded.size() <= _dim, "Inconsistent sizes.\n" )
        __ASSERT( function.coefs.size() != nb_ranks,
                  "Function's rank-coefficients have not been set.\n" )
        __ASSERT( function.basis.size() != nb_ranks,
                  "Inconsistent rank-basis size.\n" )
       
        if( _b.size() != _coefs[_dim].size() )
        {
          _b.resize( _coefs[_dim].size() );
          _A.resize( _b.size(), _b.size() );
        }
        __ASSERT( _A.size1() != _b.size(), "Wrong size for matrix _A.\n" );
        __ASSERT( _A.size2() != _b.size(), "Wrong size for matrix _A.\n" );
       
        typename T_VECTOR :: iterator i_ctarget = _b.begin();
        typename t_Expanded :: value_type 
                            :: const_iterator i_rexpI = expanded[_dim].begin();
        typename t_Expanded :: value_type
                            :: const_iterator i_rexp_end = expanded[_dim].end();
        typename t_Sizes :: value_type
                         :: const_iterator i_sizeI = sizes[_dim].begin();
        typename t_Sizes :: value_type :: const_iterator i_size_first( i_sizeI );
        size_t k1(0);
  
        // loops over ranks I
        for( size_t rI(0); i_rexpI < i_rexp_end; ++rI, ++i_rexpI, ++i_sizeI ) 
        {
          // loop over basis functions I
          typename t_Expanded :: value_type :: value_type
                              :: const_iterator i_iexp_inc = i_rexpI->begin();
          for( size_t iI(0); iI < *i_sizeI;
               ++iI, i_iexp_inc += nb_targets, ++i_ctarget, ++k1 )
          {
            __ASSERT( i_ctarget == _b.end(),        "Iterator out of range.\n" )
            __ASSERT( i_iexp_inc == i_rexpI->end(), "Iterator out of range.\n" )
            *i_ctarget = typename T_VECTOR :: value_type(0);
       
            // loops over ranks I
            typename t_Expanded :: value_type 
                                :: const_iterator i_rexpII = expanded[_dim].begin();
            typename t_Sizes :: value_type 
                             :: const_iterator i_sizeII( i_size_first );
            for( size_t rII(0), k2(0); i_rexpII < i_rexp_end;
                 ++rII, ++i_rexpII, ++i_sizeII ) 
            {
              __ASSERT( i_sizeII == sizes[_dim].end(), "Iterator out of range.\n" )
       
              // loop over basis functions II
              typename t_Expanded :: value_type :: value_type
                                  :: const_iterator i_iexpII = i_rexpII->begin();
              for(size_t iII(0); iII < *i_sizeII; ++iII, ++k2) //, ++i_A)
              {
                // Initializes element to 0 before summing over structural targets.
                _A( k1, k2 ) = 0;
       
                // loop over target values
                typename t_Expanded :: value_type :: value_type
                                    :: const_iterator i_iexpI( i_iexp_inc );
                typename T_VECTOR :: const_iterator i_starget = _targets.begin();
                typename T_VECTOR :: const_iterator i_starget_end = _targets.end();
                typename t_Weights :: const_iterator i_weight = weights.begin();
                typename t_eWeights :: const_iterator i_eweights = eweights.begin();
                for( types::t_unsigned o(0);
                     i_starget != i_starget_end;
                     ++i_starget, ++i_weight, ++i_eweights )
                {
                  __ASSERT( i_weight == weights.end(), "Iterator out of range.\n" )
                  __ASSERT( i_eweights == eweights.end(), "Iterator out of range.\n" )
       
                  // Computes first factor.
                  t_eWeights :: value_type :: const_iterator
                    i_eweight = i_eweights->begin();
                  t_eWeights :: value_type :: const_iterator
                    i_eweight_end = i_eweights->end();
                  t_Type ak1(0), ak2(0);
                  typename t_Factors :: value_type :: const_iterator i_fac;
                  for(; i_eweight != i_eweight_end;
                      ++i_eweight, ++o, ++i_iexpI, ++i_iexpII )
                     
                  {
                    __ASSERT( i_iexpI == i_rexpI->end(), "Iterator out of range.\n" )
                    __ASSERT( i_iexpII == i_rexpII->end(), "Iterator out of range.\n" )
                    __ASSERT( factors.size() <= rI * nb_targets + o, 
                              "Unexpected array-size.\n" )
                    __ASSERT( factors.size() <= rII * nb_targets + o, 
                            "Unexpected array-size.\n" )
                    __ASSERT( factors[rI * nb_targets + o].size() != D,
                              "Iterator out of range.\n" )
                    __ASSERT( factors[rII * nb_targets + o].size() != D,
                              "Iterator out of range.\n" )
                    __ASSERT( function.coefs.size() <= rI,
                              "Rank index out of range.\n" )
                    __ASSERT( function.coefs.size() <= rII,
                              "Rank index out of range.\n" )

                    if( Fuzzy::neq( *i_iexpI, 0e0 ) )
                    {
                      t_Type UI( function[ rI ] );
                      i_fac = factors[ rI * nb_targets + o].begin();
                      for( types::t_unsigned d(0); d < D; ++i_fac, ++d )
                        if( d != _dim )  UI *= (*i_fac);
                     
                      // first half of A(k,k')
                      ak1 += (*i_iexpI) * UI * (*i_eweight);
                     
                      // Computes collapsed target element.
                      if( i_sizeII == i_size_first and iII == 0 ) 
                      {
                        __ASSERT( i_ctarget == _b.end(), "Iterator out of range.\n" )
                        *i_ctarget += (*i_iexpI) * UI * (*i_starget) * (*i_eweight);
                      }
                    }
                    
                    if( Fuzzy::eq( *i_iexpII, 0e0 ) ) continue;

                    // Computes second factor.
                    t_Type UII( function[ rII ] ); 
                    i_fac = factors[ rII * nb_targets + o].begin();
                    for(types::t_unsigned d(0); d < D; ++i_fac, ++d )
                      if( d != _dim ) UII *= (*i_fac);
                    
                    // second half of A(k,k')
                    ak2 += UII * (*i_iexpII ) * (*i_eweight);
                   
                  } // end of loop over equivalent structures.
       
                  // Computes matrix element for all equiv configurations.
                  _A(k1,k2) += ak1 * ak2 * (*i_weight);
       
                } // loop over structural target values.
       
              } // loop over basis functions II
       
            } // end of loop over ranks II

          } // loop over basis functions I
          __ASSERT( i_iexp_inc != i_rexpI->end(),
                    "Unexpected iterator position.\n" )
       
        } // end of loop over ranks I
       

#       ifdef __DOHALFHALF__
          // adds regularization.
          if( Fuzzy::leq( regular_factor, 0e0 ) ) return;
          //  -- loop over ranks
          __ASSERT( _coefs[_dim].size() % 2 != 0, 
                    "Odd number of coefficients: " << _coefs[_dim].size() << ".\n" )
          typename t_Function :: t_Coefs  
                              :: const_iterator i_rcoef = function.coefs.begin();
          for( types::t_unsigned k(0); k < _coefs[_dim].size(); k+=2, ++i_rcoef )
             _A( k, k ) += regular_factor * (*i_rcoef) * (*i_rcoef) ; 
#       endif

        __DEBUGTRYEND(, "Error while creating collapsed matrix and target.\n" )
     
      } // end of functor.


  template< class T_FUNCTION > template< class T_VECTORS, class T_VECTOR > 
    typename T_FUNCTION :: t_Return
      EquivCollapse<T_FUNCTION> :: evaluate( const T_VECTORS &_coefs,
                                             const T_VECTOR &_targets ) 
      {
        namespace bl = boost::lambda;
        __DEBUGTRYBEGIN

        typename t_Factors :: const_iterator i_facs = factors.begin();
        std::vector< t_Type > values( _targets.size(), 0 );
        for( size_t r(0); r < nb_ranks; ++r )
        {
          typename std::vector< t_Type > :: iterator i_val = values.begin();
          t_eWeights :: const_iterator i_eweights = eweights.begin();
          t_eWeights :: const_iterator i_eweights_end = eweights.end();
          for(; i_eweights != i_eweights_end; ++i_val, ++i_eweights )
          {
            t_eWeights :: value_type :: const_iterator i_eweight = i_eweights->begin();
            t_eWeights :: value_type :: const_iterator i_eweight_end(i_eweights->end());
            for(; i_eweight != i_eweight_end; ++i_facs, ++i_eweight )
              *i_val += std::accumulate
                        (
                          i_facs->begin(), i_facs->end(), 1e0,
                          std::multiplies<t_Type>()
                        ) * (*i_eweight) * function[r];
          } // loop over equivalent structures.
        } // end of loop over ranks.

        t_Type result(0), dummy(0);
        opt::concurrent_loop
        (
          values.begin(), values.end(), _targets.begin(), weights.begin(),
          bl::var( result ) += ( bl::var(dummy) = bl::_2 - bl::_1,
                                 bl::var(dummy) * bl::var(dummy) * bl::_3 )
        );
#       ifdef __DOHALFHALF__
          if( Fuzzy::leq( regular_factor, 0e0 ) )
            return result / (types::t_real) _targets.size();
          typename T_VECTORS :: const_iterator i_coefs = _coefs.begin();
          typename T_VECTORS :: const_iterator i_coefs_end = _coefs.end();
          for(; i_coefs != i_coefs_end; ++i_coefs )
          {
            __ASSERT( i_coefs->size() % 2 != 0, "Odd number of coefficients.\n" )
            typename T_VECTORS :: value_type
                               :: const_iterator i_coef = i_coefs->begin();
            typename T_VECTORS :: value_type
                               :: const_iterator i_coef_end = i_coefs->end();
            typename t_Function :: t_Coefs  
                                :: const_iterator i_rcoef = function.coefs.begin();
            for(; i_coef != i_coef_end; i_coef+=2, ++i_rcoef )
            {
              types::t_real a = (*i_coef) * (*i_rcoef);
              result += regular_factor * a * a * 0.5e0; 
            }
          }
#       endif
        return result / (types::t_real) _targets.size();
        __DEBUGTRYEND(,"Error while assessing sum of squares.\n" )
      }

}