#include "NCSANSDiluteIsotropicScatter.hh"
#include "NCrystal/internal/NCMath.hh"

NCP::SANSDiluteIsotropicScatter::SANSDiluteIsotropicScatter( NC::SANSScaleFactor scalefactor,
                                                             SANSDiluteIsotropicModel themodel ) 
 : m_xsscale( scalefactor.get() * themodel.form_volume ),
    m_iofq([&themodel]() -> NC::IofQHelper {
      const int nPoints = 100000; //review this value
      const double q_min = 1e-6;
      const double Emax = 100.0;//eV
      NC::VectD Q = NC::geomspace(q_min,2.0*std::sqrt(NC::ekin2ksq(Emax)),nPoints);
      NC::VectD IofQ;
      IofQ.reserve(Q.size());
      nc_assert_always( themodel.pofq != nullptr );
      for( auto i : NC::ncrange(Q.size()))
        IofQ.push_back( themodel.pofq( Q[i] ) );
      //TODO: Move to:
      //template<class TVector, class Func>
      //inline TVector vectorTrf(const TVector&, const Func&);//
      return { Q, IofQ };
    }())
{
  nc_assert_always( themodel.form_volume > 0.0 );
}

NC::CrossSect NCP::SANSDiluteIsotropicScatter::crossSectionIsotropic(NC::CachePtr& cacheptr, NC::NeutronEnergy ekin) const
{
  double ksq = ekin.ksq();
  //Guard against zero-division (adding a crude brute-force cutoff, likely good enough).
  constexpr double ksq_tiny = 1e-250;
  if ( !( ksq > ksq_tiny ) )
    return this->crossSectionIsotropic( cacheptr, NC::NeutronEnergy{ NC::ksq2ekin(ksq_tiny) } );

  //TODO: If 2*k > largest q in m_iofq, consider extrapolation.
  return NC::CrossSect{ m_iofq.calcQIofQIntegral( ekin ) * m_xsscale / ksq };
}

NC::ScatterOutcomeIsotropic NCP::SANSDiluteIsotropicScatter::sampleScatterIsotropic( NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const
{
  double k = std::sqrt( ekin.ksq() );//NB: .k() did not work in ncrystal v3.6.0
  constexpr double k_tiny = 1e-125;
  if ( !(k>k_tiny) ) {
    ekin = NC::NeutronEnergy{ NC::ksq2ekin(k_tiny*k_tiny) };
    k = std::sqrt( ekin.ksq() );
    nc_assert_always(k>0.0);
  }
  double q = m_iofq.sampleQValue( rng, ekin );
  //FIXME: USE THE FOLLOWING:
  //  NC::CosineScatAngle mu{ NC::ncclamp(std::cos( 2.0 * std::asin( NC::ncclamp(q / ( 2.0 * k ),0.0,1.0) ) ),-1.0,1.0) };
    NC::CosineScatAngle mu{ NC::ncclamp(std::cos( 2.0 * std::asin( q / ( 2.0 * k ) ) ),-1.0,1.0) };
  return { ekin, mu };
}
