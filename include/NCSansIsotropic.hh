#ifndef NCPlugin_SansIsotropic_hh
#define NCPlugin_SansIsotropic_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

#include "NCLookUpTable.hh"
#include "NCrystal/internal/NCPointwiseDist.hh"

namespace NCPluginNamespace {

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class SansIsotropic final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static SansIsotropic createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors

    //The dummy model we are implementing is completely nonsense from a physics
    //point of view, and provides a constant sigma for wavelengths below a
    //certain cutoff value. Scatterings are isotropic and elastic.

    //Constructor gets constant cross section value, and the neutron wavelength
    //cutoff:
    SansIsotropic( const std::vector<double>& x, const std::vector<double>& f );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RandomBase& rng, double neutron_ekin ) const;

  private:
    //Data members:
    std::unique_ptr<LookUpTable> m_xs;
    std::unique_ptr<NC::PointwiseDist> m_sansQDist;
  };

}
#endif
