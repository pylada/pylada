//
//  Version: $Id$
//
#ifndef _DARWIN_EVALUATOR_H_
#define _DARWIN_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <string>
#include <algorithm>

#include <eo/eoGenOp.h>

#include "opt/types.h"
#include "opt/opt_function_base.h"
#include "taboos.h"
#include "loadsave.h"

#include "gatraits.h"

namespace GA 
{
  template< class T_INDIVIDUAL>
  class Evaluator
  {
    public:
      typedef T_INDIVIDUAL   t_Individual;
    protected:
      typedef typename t_Individual :: t_IndivTraits        t_IndivTraits;
      typedef typename t_IndivTraits :: t_Object            t_Object;
      typedef typename t_IndivTraits :: t_QuantityTraits    t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity       t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_IndivTraits :: t_VA_Traits         t_VATraits;
      typedef typename t_VATraits :: t_QuantityGradients    t_QuantityGradients;
      typedef typename t_VATraits :: t_Type                 t_VA_Type;

    protected:
      t_Individual *current_individual;
      t_Object *current_object;

    public:
      Evaluator() {};
      ~Evaluator() {}

    public:
      // Should load t_Individual and funtional related stuff
      // except attributes in <GA > tag
      bool Load ( std::string const &_f )
      {
        TiXmlDocument doc( _f.c_str() ); 
        TiXmlHandle docHandle( &doc ); 
        if  ( !doc.LoadFile() )
        { 
          std::cerr << doc.ErrorDesc() << std::endl; 
          throw "Could not load input file in CE::Evaluator ";
        } 

        TiXmlElement *child = docHandle.FirstChild("Job").Element();
        if (not child)
          return false;
        return Load(*child);
      }
      bool Load ( const TiXmlElement &_node ) { return true; }
      // Load and Save individuals
      // _type can be either GA::LOADSAVE_SHORT or
      // GA::LOADSAVE_LONG. Results are save as latter, and GA
      // internal stuff as the former. You need both only if the
      // ga object and the user-expected result object are different (say
      // bitstring versus a decorated lattice structure )
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type ) {return true;};
      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const {return true;};
      // attributes from <GA > tag in input.xml are passed to this
      // function from Darwin::Load(...)
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      // returns a pointer to an eoOp object
      // pointer is owned by GA::Darwin::eostates !!
      // don't deallocate yourself
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el ) { return NULL; };
      // returns a pointer to a eoF<bool> object
      // pointer is owned by GA::Darwin::eostates !!
      // don't deallocate yourself
      eoF<bool>* LoadContinue(const TiXmlElement &_el ) { return NULL; }
      // returns a pointer to a eoMonOp<const t_Individual> object
      // pointer is owned by GA::Darwin::eostates !!
      // don't deallocate yourself
      Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el ) { return NULL; }
      // Initializes object before call to functional 
      // i.e. transforms t_Individual format to
      // function::Function<...>::variables format if necessary
      bool initialize( t_Individual &_indiv ) {return false; };
      // Called before objective function is evaluated
      // must return a void pointer to functional
      void init( t_Individual &_indiv )
      {
        current_individual = &_indiv;
        current_object     = &_indiv.Object();
      };
      void evaluate() {};
      // Override the next three functions only if VA Minimization is implemented
      void evaluate_gradient( t_QuantityGradients& _grad )
        { Traits::zero_out( _grad ); }
      void evaluate_with_gradient( t_QuantityGradients& _grad )
      {
        evaluate_gradient( _grad );
        evaluate();
      }
      void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos) 
        { Traits :: zero_out( _grad[_pos] ); }
      //! \brief Submits individuals to history, etc, prior to starting %GA
      //! \details initializes the endopoints of a convex-hull, for instance.
      //! Presubmitted individuals are not put into the population.
      //! \note GA::Evaluator does not know about Traits::GA (or affiliated).
      //! Hence it does not know about Traits::GA::t_Population. So we input this
      //! population an eoPop<t_Individual> directly.
      void presubmit( std::list<t_Individual>& ) { return; }
  };

}

#endif // _EVALUATOR_H_
