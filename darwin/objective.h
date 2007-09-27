//
//  Version: $Id$
//
#ifndef _MULTIOB_OBJECTIVE_H_
#define _MULTIOB_OBJECTIVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/convex_hull.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif
#include <print/xmg.h>
#include <print/manip.h>

#include "gatraits.h"
#include "loadsave.h"

namespace Objective
{

  template< class T_GA_TRAITS,
            class T_QUANTITY_TRAITS = typename T_GA_TRAITS :: t_QuantityTraits,
            class T_VA_TRAITS = typename T_GA_TRAITS :: t_VA_Traits >
  class Base
  {
    public: 
      typedef T_GA_TRAITS t_GATraits;
      typedef T_QUANTITY_TRAITS t_QuantityTraits;
      typedef T_VA_TRAITS t_VA_Traits;
    protected:
      typedef typename t_GATraits :: t_Evaluator           t_Evaluator;
      typedef typename t_GATraits :: t_Individual          t_Individual;
      typedef typename t_GATraits :: t_Fitness             t_RealFitness;
      typedef typename Fitness::GetScalar< t_RealFitness, 
                                           t_QuantityTraits::is_scalar > :: t_result t_Fitness;
      typedef typename t_QuantityTraits :: t_Quantity       t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_VA_Traits :: t_QuantityGradients   t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type                t_VA_Type;
      typedef GA::SaveObject<t_GATraits>                    t_SaveOp;
      typedef GA::LoadObject<t_GATraits>                    t_LoadOp;
    protected:
      static const t_Individual *current_indiv;
      static t_Fitness fitness;
    public:
      bool const static is_scalar = t_QuantityTraits :: is_scalar;
      bool const static is_vector = t_QuantityTraits :: is_vector;
    public:
      Base() {}
      virtual ~Base() {}

      void static init(const t_Individual& _indiv);
      virtual const t_Fitness& operator()( const t_Quantity& ) = 0;
      virtual void evaluate_gradient( const t_Quantity &,
                                      t_QuantityGradients &,
                                      t_ScalarQuantity *) = 0;
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *) = 0;
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients&,
                                               types::t_unsigned) = 0;
      virtual bool is_valid() const = 0;
      virtual std::string what_is() const = 0;
      virtual std::string print() const = 0;
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op) { return true; };
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) { return true; };
      virtual bool does_store() const { return false; }
  };
  template< class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
  void Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: init( const t_Individual & _indiv)
    { Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = &_indiv; }
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    const typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Individual* 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = NULL;
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Fitness 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: fitness;
 
  template< class T_TYPE, bool IS_VECTORIAL > struct fork;

  template <class T_GA_TRAITS >
    struct Types
    {
      template< class T_TYPE, bool is_vectorial > friend struct fork;
      public:
        typedef T_GA_TRAITS t_GATraits;
      protected:
        typedef typename t_GATraits::t_Evaluator            t_Evaluator;
        typedef typename t_GATraits::t_Individual           t_Individual;
        typedef typename t_GATraits::t_QuantityTraits       t_QuantityTraits;
        typedef typename t_QuantityTraits::t_ScalarQuantity t_ScalarQuantity;

      public:
        typedef Base< t_GATraits, 
                      typename Traits::Quantity< t_ScalarQuantity >,
                      typename Traits::VA< std::vector<t_ScalarQuantity>,
                                           t_ScalarQuantity > >           Scalar;
        typedef Base< t_GATraits >  Vector;



        static Vector* new_from_xml( const TiXmlElement &_node )
          { fork<Types, t_QuantityTraits::is_vector > s; return s( _node ); }
      protected:
        static Scalar* scalar_from_xml( const TiXmlElement &_node );
        static Vector* vector_from_xml( const TiXmlElement &_node );
    };

  template< class T_GA_TRAITS >
  class Maximize : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Evaluator         t_Evaluator;
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;

    protected:
      using t_Base :: fitness;

    public:
      Maximize() {}
      Maximize( const Maximize &) {}
      virtual ~Maximize() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val)
       { fitness = -_val; return fitness; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result -= *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result -= *i_grad;
        return -_val;
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return -_grad[_n]; }
      bool is_valid() const { return true; }
      virtual std::string what_is() const { return " Maximize"; }
      virtual std::string print() const { return ""; }
  };
  template< class T_GA_TRAITS >
  class Minimize : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;

    protected:
      using t_Base :: fitness;

    public:
      Minimize() {}
      Minimize( const Minimize &) {}
      virtual ~Minimize() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val)
         { fitness = _val; return fitness; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result += *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result += *i_grad;
        return _val;
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return _grad[_n]; }
      bool is_valid() const { return true; }
      virtual std::string what_is() const { return " Minimize"; }
      virtual std::string print() const { return ""; }
  };
  template< class T_GA_TRAITS >
  class Target : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_GATraits :: t_IndivTraits       t_IndivTraits;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;

    protected:
      using t_Base :: fitness;

    protected:
      t_ScalarQuantity target; 
    public:
      Target( t_ScalarQuantity _target ) : target( _target ) {}
      Target( const Target &_c ) : target( _c.target ) {}
      virtual ~Target() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val)
        { fitness = std::abs( _val - target ); return fitness; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
        return std::abs(target - _val);
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) 
      {
        return ( _val > target ) ? _grad[_n]: -_grad[_n];
      }
      bool is_valid() const { return true; }
      virtual std::string print() const { return ""; }
      virtual std::string what_is() const
      { 
        std::ostringstream sstr;
        sstr << " Target (" << target << ")";
        return sstr.str();
      }
  };
  template< class T_GA_TRAITS >
  class ConvexHull : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef opt::ConvexHull::Base<t_Individual>         t_ConvexHull;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;

    protected:
      t_ConvexHull convexhull; 
      mutable bool valid;
      using t_Base :: fitness;
      using t_Base :: current_indiv;

    public:
      ConvexHull() : valid(true) {}
      ConvexHull( const ConvexHull &_c ) : convexhull( _c.convexhull ), valid(true) {}
      virtual ~ConvexHull() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val);
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
        { evaluate_with_gradient( _val, _grad, _i_grad ); }

      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *);
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      bool is_valid() const
      {
        if ( valid ) return true;
        valid = true;
        return false;
      }
      virtual std::string what_is() const { return " Convex-Hull"; }
      virtual std::string print() const { return convexhull.print(); }
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
        { return convexhull.Save( _node, _op ); };
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op)
        { return convexhull.Load( _node, _op ); };
      virtual bool does_store() const { return true; }
  };
  template< class T_GA_TRAITS >
  const typename ConvexHull<T_GA_TRAITS> :: t_Fitness&
    ConvexHull<T_GA_TRAITS> :: operator()(const t_Quantity& _val)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
    
      if ( _val >= base )
       { fitness = _val - base;  return fitness; }
    
      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      fitness = 0.0;
      return fitness;
    }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_ScalarQuantity
    ConvexHull<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients &_grad,
                                                       t_VA_Type *_i_grad)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
      typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
      t_VA_Type *i_grad_result = _i_grad;
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result += *i_grad - gradient;
      
      if ( _val >= base ) return _val - base;

      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      return 0.0;
    }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_VA_Type
    ConvexHull<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                      t_QuantityGradients &_grad,
                                                      types::t_unsigned _n) 
    {
      t_Quantity x = current_indiv->get_concentration();
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      return _grad[_n] - gradient;
    }


  template<class T_GA_TRAITS >
  class Container : public Types<T_GA_TRAITS> :: Vector
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef Types<t_GATraits>                          t_ObjectiveType;
      typedef typename t_ObjectiveType :: Vector          t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_GATraits :: t_IndivTraits       t_IndivTraits;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_ObjectiveType :: Scalar          t_Objective;
      typedef std::vector< t_Objective* >                 t_Objectives;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;

    protected:
      using t_Base :: fitness;
      t_Objectives objectives;

    public:
      Container() {}
      Container( const Container &_c ) : objectives(_c.objectives) {}
      virtual ~Container () 
      {
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.end();
        for(; i_objective != i_end; ++i_objective ) delete *i_objective;
        objectives.clear();
      }

      bool is_valid() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->is_valid() ) return false;
        return true;
      }
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->Save( _node, _op ) ) return false;
        return true;
      }
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) 
      {
        typename t_Objectives::iterator i_objective = objectives.begin();
        typename t_Objectives::iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->Restart( _node, _op ) ) return false;
        return true;
      }
      virtual bool does_store() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( (*i_objective)->does_store() ) return true;
        return false;
      }
      virtual std::string print() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        std::ostringstream sstr;
        for(; i_objective != i_end; ++i_objective )
          sstr << (*i_objective)->print();
        return sstr.str();
      }
  };

  template<class T_GA_TRAITS >
  class LinearSum : public Container<T_GA_TRAITS>
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef Types<t_GATraits>                      t_ObjectiveType;
      typedef Container<t_GATraits>                  t_Base;
      typedef typename t_Base :: t_Fitness           t_Fitness;
      typedef typename t_Base :: t_Individual        t_Individual;
      typedef typename t_Base :: t_QuantityTraits    t_QuantityTraits;
      typedef typename t_Base :: t_Quantity          t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity    t_ScalarQuantity;
      typedef typename t_Base :: t_Objective         t_Objective;
      typedef std::vector< t_Objective* >            t_Objectives;
      typedef typename t_Base :: t_VA_Traits         t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>             t_SaveOp;
      typedef GA::LoadObject<t_GATraits>             t_LoadOp;

    protected:
      using t_Base :: fitness;
      using t_Base :: objectives;
      std::vector< t_ScalarQuantity > coefs;

    public:
      LinearSum() {}
      LinearSum ( const LinearSum &_c ) : t_Base(_c), coefs(_c.coefs) {}
      virtual ~LinearSum() {}

      void add( t_Objective *_objective, t_ScalarQuantity _coef )
      {
        if ( not _objective ) return;
        coefs.push_back( _coef );
        objectives.push_back( _objective );
      }
      virtual const t_Fitness& operator()(const t_Quantity& _val);
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *);
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      virtual std::string what_is() const
      {
        std::ostringstream sstr;
        sstr << "LinearSum begin{ ";
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        for(; i_objective != i_end; ++i_objective, ++i_coef )
          sstr << (*i_objective)->what_is() << "[" << *i_coef << "] ";
        sstr << "} end"; 
        return  sstr.str();
      }
     
  };
  
  template< class T_GA_TRAITS >
    const typename LinearSum<T_GA_TRAITS>::t_Fitness&
      LinearSum<T_GA_TRAITS> :: operator()( const t_Quantity& _val ) 
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );

        t_ScalarQuantity inter = 0;
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();
        fitness.clear();
        for(; i_objective != i_end and i_val != i_val_end; ++i_objective, ++i_coef, ++i_val )
        {
          double r = (*i_objective)->operator()( *i_val );
          fitness.push_back( r );
          inter += ( *i_coef ) * r;
        }
          
        fitness = inter;
        return fitness;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_ScalarQuantity
      LinearSum<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                                    t_QuantityGradients &_grad,
                                                                    t_VA_Type *_i_grad)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        t_ScalarQuantity results = 0.0;
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type *i2 = _i_grad + _val.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _val.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _val.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          results +=   (*i_objective)->evaluate_with_gradient( *i_val, *i_grad, const_grad_result )
                     * (*i_coef);
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
          
        return results;
      };
  template< class T_GA_TRAITS >
    void
      LinearSum<T_GA_TRAITS> ::  evaluate_gradient( const t_Quantity &_val,
                                                                t_QuantityGradients &_grad,
                                                                t_VA_Type *_i_grad)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type *i2 = _i_grad + _val.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _val.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _val.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          (*i_objective)->evaluate_gradient( *i_val, *i_grad, const_grad_result );
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_VA_Type
      LinearSum<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity & _val,
                                                                   t_QuantityGradients& _grad,
                                                                   types::t_unsigned _n)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type result = t_VA_Type(0);
        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
          result +=   ( *i_coef ) 
                    * (*i_objective)->evaluate_one_gradient( *i_val, *i_grad, _n );

        return result;
      };

  template< class T_GA_TRAITS >
     typename Types<T_GA_TRAITS> :: Scalar*
      Types<T_GA_TRAITS> :: scalar_from_xml( const TiXmlElement &_node )
      {
        if ( not &_node ) return NULL;
        std::string str = "minimize"; 
        std::string name = Print::lowercase(_node.Value());
        if (    name.compare("objective") == 0 
             or name.compare("method") == 0 )
        {
          if ( _node.Attribute( "type" ) )
            str = Print::lowercase(_node.Attribute( "type" ));
        }
        else if ( _node.Attribute("objective") )
          str = _node.Attribute( "objective" );
        if ( str.compare("convexhull") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: ConvexHull" << Print::endl;
          return new ConvexHull<t_GATraits>;
        }
        else if ( str.compare("minimize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Minimize" << Print::endl;
          return new Minimize<t_GATraits>;
        }
        else if ( str.compare("maximize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Maximize" << Print::endl;
          return new Maximize<t_GATraits>;
        }
        else if (str.compare("target") == 0 )
        {
          if( _node.Attribute("target") )
          {
            double d; _node.Attribute("target", &d );
            Print::xmg << Print::Xmg::comment
                       << "Objective: Target (" << d << ")" << Print::endl;
            return new Target<t_GATraits>( (types::t_real) d );
          }
        }
        if ( _node.FirstChildElement( "Objective" ) )
         return scalar_from_xml( *_node.FirstChildElement( "Objective" ) ); 

        return NULL;
      }
    template< class T_GA_TRAITS >
       typename Types<T_GA_TRAITS> :: Vector*
        Types<T_GA_TRAITS> :: vector_from_xml( const TiXmlElement &_node )
        {
          if ( not &_node ) return NULL;
          std::string str = "minimize"; 
          std::string name = Print::lowercase(_node.Value());
          if (    name.compare("objective") == 0 
               or name.compare("method") == 0 )
          {
            if ( _node.Attribute( "type" ) )
              str = Print::lowercase(_node.Attribute( "type" ));
          }
          if ( Vector::t_QuantityTraits::is_vector ) // and str.compare("LinearSum") == 0 )
          {
            LinearSum<T_GA_TRAITS> *linear = new LinearSum<T_GA_TRAITS>;
            if ( not linear ) 
            {
              std::cerr << "Mememory Pb when creating LinearSum multi-objective" << std::endl;
              return NULL;
            }
            const TiXmlElement *child = child->FirstChildElement("Objective");
            for(; child; child = child->NextSiblingElement("Objective") )
            {
              Scalar* scalar = scalar_from_xml( *child );
              if ( not scalar ) continue;
              double d = 0.0;
              if ( not child->Attribute("coef", &d) ) d = 1.0;
              linear->add( scalar, t_ScalarQuantity(d) );
            }
            return linear;
          }

          if ( _node.FirstChildElement( "Objective" ) )
           return new_from_xml( *_node.FirstChildElement( "Objective" ) ); 

          return NULL;
        }

   
    // bullshit class since there is no such thing as partial specialization of functions!!
    template< class T_TYPE >
      struct fork<T_TYPE, true>
      {
        typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
        {
          typename T_TYPE :: Vector* result;
          result = T_TYPE::vector_from_xml( _node );
          if( not result ) std::cerr << "Could not create multi-objective..." << std::endl;
          return result;
        }
      };
    template< class T_TYPE > 
      struct fork<T_TYPE, false>
      {
        typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
        {
          return T_TYPE::scalar_from_xml( _node );
        }
      };

  }  // namespace Objective
#endif //  _MULTIOB_OBJECTIVE_H_