//
//  Version: $Id$
//
#ifndef _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_ 
# define _PESCAN_INTERFACE_DIPOLE_ELEMENTS_H_
# ifdef _DIRECTIAGA 

# ifdef HAVE_CONFIG_H
#   include <config.h>
# endif

# include<vector>
# include<complex>

# include<opt/types.h>

namespace LaDa
{
  namespace Crystal 
  {
    class Structure;
  }

  namespace Pescan
  {
    class BandGap;

    //! Holds dipole moments.
    struct Dipole
    {
      friend class boost::serialization::access;
      typedef std::pair<size_t, size_t> t_Band2Band;
      //! The complex dipole moment.
      std::complex<types::t_real> r[3];
      //! Indices of bands in transition.
      t_Band2Band band2band;
      private:
        //! Serializes a dipole moment.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          {  _ar & r; _ar & band2band; }
    };
    //! Prints a dipole moment.
    std::ostream& operator<<( std::ostream &_stream, const Dipole& );
    //! \brief Computes the transition dipole elements between valence and conduction bands.
    //! \details The transition dipole element is given by the trace of
    //!           \f$-i[\hat{\matcal{H}}, \hat{r}]\f$. See <A
    //!           HREF="http://dx.doi.org/10.1103/PhysRevA.3.1242"> Anthony
    //!           F. Starace, PRA <STRONG>3</STRONG>, 1242 (1977). </A> for
    //!           details of the derivation. 
    //! \todo The current implementation expects the system to be spin
    //!       degenerate, e.g. Kramer's degeneracy.
    void dipole_elements( std::vector< Dipole > &_dipoles,
                          const BandGap& _bandgap,
                          const types::t_real _degeneracy = types::tolerance );
    //! Returns the valence-conduction the norm of the band dipole elements.
    types::t_real oscillator_strength( const std::vector<Dipole> &_dipoles );
    //! Returns the valence-conduction the norm of the band dipole elements.
    types::t_real oscillator_strength( const BandGap& _bandgap,
                                       const types::t_real _degeneracy = types::tolerance,
                                       bool _print = false );

  }
} // namespace LaDa

# endif
#endif
