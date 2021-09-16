#ifndef NCPlugin_SansIsotropic_hh
#define NCPlugin_SansIsotropic_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

#include "NCLookUpTable.hh"
#include "NCrystal/internal/NCPointwiseDist.hh"

namespace NCPluginNamespace {

  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class SansIsotropic final : public NC::MoveOnly {
  public:
    SansIsotropic( const std::vector<double>& Q, const std::vector<double>& I );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    std::unique_ptr<LookUpTable> m_xs;
    std::unique_ptr<NC::PointwiseDist> m_sansQDist;
  };

}
#endif
