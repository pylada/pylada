//
//  Version: $Id$
//
#ifndef _RANKING_H_
#define _RANKING_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

/** \ingroup Genetic
 * @{ */
//! \brief Scales the fitness of an individual according to the whole population
//! \details There are two kinds of scaling: one were there is some ranking
//!    done, one were fitnesses are merely scaled. Basically, in this namespace are all
//!    population dependent fitness operations. This may mean a Pareto Ranking,
//!    but it could also be niching. These operations are put together mostly
//!    because they take place at the same place, eg once all individuals have
//!    been evaluted and their objective %functions computed. 
//! \note For multi-objective fitnesses (eg vectorial), it may be tricky to
//!       know where the results of a scaling operation is stored (in the
//!       vectorial or scalar fitness?) and which it uses to do the scaling.
//!       These technicalities are duly noted for each scaling functor.
namespace Scaling
{

  //! \brief Base class for population dependent fitness operations
  //! \details A pure virtual operator(t_Population&) is defined by the base
  //!   class. 
  template<class T_GATRAITS> 
    class Base : public eoUF< typename T_GATRAITS::t_Population&, void> 
  {
    public:
    //! Should return a string defining the scaling operation
    virtual std::string what_is() const = 0;
  };

  //! Returns a scaling from an XML input
  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_from_xml( const TiXmlElement &_node );

  //! Returns a niching from an XML input
  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node );

  //! \brief A scaling which calls upon other to do the job
  //! \details Merely implements general subroutines around a container of
  //!          Scaling::Base pointers. Note that since scaling operators
  //!          generally store their result in the individuals fitness. A
  //!          container implies that scaling n+1 is applied upon the result of
  //!          scaling n. Note however that this particular Scaling object does not
  //!          even know about Fitnesses.
  //! \warning The pointers are owned by this class. They will be destroyed by
  //!          this class
  //! \note In the case of (vectorial) multi-objective %GA, it is expected that
  //!       the end product is \a scalar fitness. As such, whatever the last
  //!       ranking functor in this container may act upon, it must store its
  //!       result as a scalar fitness. Otherwise the whole operation is vain.
  template<class T_GATRAITS>
  class Container : public Base<T_GATRAITS>
  {
    public:
        typedef T_GATRAITS t_GATraits; //!< Contain all %GA types
      protected:
        //! Type of the population
        typedef typename t_GATraits :: t_Population t_Population;
        //! Type of the base class
        typedef Base<t_GATraits> t_Base;
        //! Type of the containe of scaling
        typedef std::list<t_Base*> t_Container;

      protected:
        t_Container rankers; //!< A container of Scaling pointers

      public:
        //! Constructor
        Container() {} 
        //! Copy Constructor
        Container( const Container &_c ) : rankers(_c.rankers) {}
        //! Destructor
        ~Container();

        //! Puts \a _ptr at end of container
        void push_back( t_Base* _ptr) { if( _ptr ) rankers.push_back( _ptr ); }
        //! Puts \a _ptr at beginning of container
        void push_front( t_Base* _ptr) { if( _ptr ) rankers.push_front( _ptr ); }

        //! \brief Removes the front scaling and returns a pointer to it
        //! \details It is now your job to delete the returned pointer
        Base<t_GATraits>* pop_front(); 
        //! Applies each Scaling object in turn
        void operator()( t_Population & _pop );

        //! Returns the number of stored Scaling pointers
        types::t_unsigned size() const { return rankers.size(); }

        //! Returns a string defining each stored Scaling object
        virtual std::string what_is() const;
    };

    /** \brief Defines a Niching operation \f$\mathcal{F}(\sigma)
               = \mathcal{F}(\sigma)/\sum_{\{\sigma_j\}_{(n)}}  \mathcal{D}(\sigma,
               \sigma_j)\f$.
        \details Niching is an operation which attempts to reduce the fitness of
                 individuals which are, somehow to close to one another. Using
                 a sharing %function \f$\mathcal{S}(\sigma, \sigma_j)\f$ which
                 judges this closeness, it implements a population dependent
                 fitness.
                 \f[
                     \mathcal{F}(\sigma) = 
                         \frac{\mathcal{F}(\sigma)}
                              {\sum_{\{\sigma_j\}} \mathcal{S}(\sigma, \sigma_j)}
                  \f]
                  Where \f$\mathcal{F}(\sigma)\f$ is the fitness of individual
                  \f$\sigma\f$, \f$\{\sigma\}\f$ represents the current
                  population. The sharign %function can be pretty much anything
                  you want, but should reflect how genetically similar two
                  individuals are.
        \note Whether for single or multi-objective %GA, this scaling functor
              acts upon the \b scalar fitness only. 
        \note \f$\mathcal{S}(\sigma_i, \sigma_i)\f$ automatically defaults to 1.
  */
  template<class T_SHARING>
  class Niching : public Base<typename T_SHARING::t_GATraits >
  {
    public:
      typedef T_SHARING t_Sharing; //!< The Sharing %function
      typedef typename t_Sharing :: t_GATraits t_GATraits; //!< All %GA %types
    protected:
      //! The type of the base class
      typedef Base< t_GATraits > t_Base;
      //! The type of the individuals
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! The type of the fitness
      typedef typename t_GATraits :: t_Fitness t_Fitness;
      //! The type of the population
      typedef typename t_GATraits :: t_Population t_Population;
      //! The type of the quantity of the fitness
      typedef typename t_Fitness  :: t_Quantity t_FitnessQuantity;
      //! The type of the scalar fitness
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      //! The type of scalar quantity of the fitness
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      //! An internal type which holds the scaling of each individual
      typedef std::vector< t_ScalarFitnessQuantity >  t_Column;

    protected:
      t_Sharing sharing; //!< The sharing %function
      t_Column sums; //!< The vector \f$\sum_j\mathcal{S}(\sigma_i, \sigma_j)\f$

    public:
      //! Constructor
      Niching() : t_Base(), sharing() {}
      //! Constructor and Initialiser 
      Niching( const TiXmlElement &_node ) : t_Base(), sharing( _node ) {}
      //! Copy Constructor
      Niching( const Niching &_n) : t_Base(_n), sharing(_n.sharing), sums(_n.sums) {}
      //! Destructor
      ~Niching() {}
  
      //! Scales down the fitness of each individuals in \a _pop according to
      //! its closeness to other individuals
      void operator()(t_Population& _pop);
      //! Loads the sharing operator from XML
      bool Load( const TiXmlElement &_node )
        { return sharing.Load(_node); }
      //! Returns a string describing the niching
      virtual std::string what_is() const 
        { return "Scaling::Niching[ " + sharing.what_is() + " ]"; }


    protected:
      //! \brief Computes the vector \f$\sum_j\mathcal{S}(\sigma_i, \sigma_j)\f$
      //! \details \f$\mathcal{S}(\sigma_i, \sigma_i)\f$ automatically defaults
      //!         to 1.
      void map(const t_Population &_pop);
  };

  /** \brief Ranks individuals according to how dominating they are
      \details Let \f$\mathcal{F}^{v}_t(\sigma)\f$ denote the component
               \f$t\in[0,N[\f$ of <I>N</I>-dimensional fitness 
               \f$\mathcal{F}^{v}(\sigma)\f$ of an individual \f$\sigma\f$.
               Also, let \f$\mathcal{F}(\sigma)\f$ be the \b scalar end-product
               fitness of an individual \f$\sigma\f$.  An individual
               \f$\sigma_i\f$ dominates another individual \f$\sigma\f$ if
               and only if 
               \f[
                   \mathcal{F}^v(\sigma_i) \succeq \mathcal{F}^v(\sigma_j)\quad
                   \Leftrightarrow\quad \forall t\in[0,N[,\ \mathcal{F}_t(\sigma_i)
                   \geq \mathcal{F}_t(\sigma_j).
               \f]
               This definition of Pareto ordering is defined by Fitness::Base
               (both for vectorial and scalar flavors). ParetoRanking sets the
               end-product fitness of each individual to
               \f[
                 \mathcal{F}(\sigma) = \sum_j\mathcal{F}^v(\sigma) \succeq
                                             \mathcal{F}^v(\sigma_j)\ ?\ -1 :\ 0 
               \f]
               The negation is there because the default is to minimize. Also,
               according to the definition, individuals dominate themselves.
      \note For a multi-objective %GA, ParetoRanking acts upon the \b vectorial
            fitness and sets the \b scalar fitness. For single-objective its
            acts, by default, upon the \b scalar fitness.
  */
  template<class T_GATRAITS>
  class ParetoRanking : public Base< T_GATRAITS >
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA %type
    protected:
      //! The type of the base class
      typedef Base< t_GATraits > t_Base;
      //! The type of the individuals
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! The type of the fitness
      typedef typename t_GATraits :: t_Fitness t_Fitness;
      //! The type of the population
      typedef typename t_GATraits :: t_Population t_Population;
      //! The type of the quantity of the fitness
      typedef typename t_Fitness  :: t_Quantity t_FitnessQuantity;
      //! The type of the scalar fitness
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      //! The type of scalar quantity of the fitness
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      //! An internal type which holds the number of individuals each
      //! individual dominates
      typedef std::vector< types::t_unsigned >       t_Column;

    protected:
      /** The vector \f$\sum_j\mathcal{F}^v(\sigma) \succeq
          \mathcal{F}^v(\sigma_j)\ ?\ -1 :\ 0  \f$ */
      t_Column sums;

    public:
  
      //! Constructor
      ParetoRanking() : t_Base() {}
      //! Constructor and Initialiser 
      ParetoRanking( const TiXmlElement &_node ) : t_Base() {}
      //! Copy Constructor
      ParetoRanking( const ParetoRanking &_n) : t_Base(_n), sums(_n.sums) {}
      //! Destructor
      ~ParetoRanking() {}
  
      /** Sets the scalar fitnes to the number of dominated individual,
          \f$ \mathcal{F}(\sigma) = \sum_j\mathcal{F}^v(\sigma) \succeq
              \mathcal{F}^v(\sigma_j)\ ?\ -1 :\ 0 \f] */
      void operator()(t_Population& _pop);

      //! Does nothing
      bool Load( const TiXmlElement &_node )  { return true; }
      //! Returns "Ranking::Pareto"
      virtual std::string what_is() const 
        { return "Ranking::Pareto"; }
     
    protected:
      /** \brief Computes the vector \f$\sum_j\mathcal{F}^v(\sigma) \succeq
          \mathcal{F}^v(\sigma_j)\ ?\ -1 :\ 0  \f$ */
      void map(const t_Population &_pop);
  };


  //! \brief Contains all sharing %function
  //! \details A sharing %function is a functor which defines a scaling from the
  //! distance between two individuals. Generally, the sharing between an
  //! individual and itself should be 1. Furthermore two individuals which are
  //! more distant than a given cutoff are said not to share (eg sharing %function
  //! is 0). Other than that, use your imagination.
  namespace Sharing
  {
    /** \brief Implements a triangular sharing %function 
     *  \details With \a T_DISTANCE a distance functor
     *  \f$|\sigma-\sigma_j|\f$, the triangular %function is simply defined as
     *  \f[ 
     *      \mathcal{S} = |\sigma_i-\sigma_j| > d_0\ ?\ 0:\ 1 -
     *                    \left(\frac{|\sigma_i-\sigma_j|}{d_0}\right)^\alpha
     *  \f]
     */
    template <class T_DISTANCE>
    class Triangular 
    {
      public:
        //! Type of the distance functor
        typedef T_DISTANCE t_Distance; 
        //! All %GA %types
        typedef typename t_Distance :: t_GATraits t_GATraits;
      protected:
        //! Type of the individuals
        typedef typename t_GATraits :: t_Individual t_Individual;
        //! Type of the fitness
        typedef typename t_GATraits :: t_Fitness t_Fitness;
        //! Type of the \b scalar fitness
        typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
        //! Type of the quantity of the \b scalar fitness
        typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

      protected:
        t_Distance distance; //!< An instance of the distance functor
        types::t_unsigned alpha; //!< The \f$\alpha\f$ parameter
        t_ScalarFitnessQuantity d_0; //!< the cutoff parameter \f$d_0\f$

      public:
        //! Constructor
        Triangular() : distance(), alpha(1), d_0(1)  {}
        //! Constructor and initialiser
        Triangular  ( const TiXmlElement &_node ) 
                   : distance(), alpha(1), d_0(1) { Load(_node); }
        //! Copy Constructor
        Triangular   (const Triangular &_t )
                   : distance(_t.distance), alpha(_t.alpha), 
                     d_0(_t.d_0)  {}
        //! Destructor
        ~Triangular() {}

        //! Loads the sharing %function from XML
        bool Load( const TiXmlElement &_node );
        /** returns \f$ \mathcal{S} = |\sigma_i-\sigma_j| > d_0\ ?\ 0:\ 1 -
                        \left(\frac{|\sigma_i-\sigma_j|}{d_0}\right)^\alpha \f$ */
        t_ScalarFitnessQuantity operator()(const t_Individual& _i1, const t_Individual& _i2 ) const;

        //! Returns a string defining the sharing %function
        std::string what_is() const;
    };  
  }

  //! \brief Defines distance binary functor bewteen individuals
  //! \details A distance is any binary functor of two individuals which returns
  //!          a \b scalar quantity. That said, it should somehow relate to
  //!          distances... In other words it should be larger for two
  //!          individuals which share a large part of their genotype or of
  //!          their phenotype (whenever the two are distincts). 
  namespace Distance
  {
    /** \brief Defines a generalized Hamming distance \f$|\sigma_i-\sigma_j|=\sum_t
       |q_{t,\sigma_i} - q_{t,\sigma_j}|\f$ */
    template<class T_GATRAITS>
    class GeneralHamming
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< All %GA %types
      protected:
        //! Type of the individuals
        typedef typename t_GATraits :: t_Individual t_Individual;
        //! Type of the object characterising the individuals
        typedef typename t_GATraits :: t_Object t_Object;
        //! Type of the \b scalar fitness
        typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
        //! Type of the  quantity of the \b scalar fitness
        typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

      public:
        //! Returns \f$|\sigma_i-\sigma_j|=\sum_t |q_{t,\sigma_i} - q_{t,\sigma_j} |\f$
        t_ScalarFitnessQuantity operator()( const t_Individual &_i1, const t_Individual &_i2) const;
        //! Returns "Distance::GeneralHamming"
        std::string what_is() const { return "Distance::GeneralHamming"; }
    }; 
    /** \brief Defines the Hamming distance \f$|\sigma_i-\sigma_j|=\sum_t
     *  (q_{t,\sigma_i} = q_{t,\sigma_j})\ ?\ 1:\ 0\f$ */
    template<class T_GATRAITS>
    class Hamming
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< All %GA %types
      protected:
        //! Type of the individuals
        typedef typename t_GATraits :: t_Individual t_Individual;
        //! Type of the object characterising the individuals
        typedef typename t_GATraits :: t_Object t_Object;
        //! Type of the \b scalar fitness
        typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
        //! Type of the  quantity of the \b scalar fitness
        typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

      public:
        //! Returns \f$\sum_t (q_{t,\sigma_i} = q_{t,\sigma_j})\ ?\ 1:\ 0 \f$
        t_ScalarFitnessQuantity operator()( const t_Individual &_i1, const t_Individual &_i2) const;
        //! Returns "Distance::Hamming"
        std::string what_is() const { return "Distance::Hamming"; }
    }; 

  } // namespace Scaling

} // namespace Distance
#include "scaling.impl.h"

/* @} */

#endif // _RANKING_H_
