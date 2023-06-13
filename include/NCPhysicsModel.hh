#ifndef NCPlugin_PhysicsModel_hh
#define NCPlugin_PhysicsModel_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/internal/NCIofQHelper.hh"
namespace NCPluginNamespace {

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class PhysicsModel final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static PhysicsModel createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors

    //The dummy model we are implementing is completely nonsense from a physics
    //point of view, and provides a constant sigma for wavelengths below a
    //certain cutoff value. Scatterings are isotropic and elastic.

    //Constructor gets constant cross section value, and the neutron wavelength
    //cutoff:
    PhysicsModel( double radius, double sld, double sld_solvent);

    // provides IQ for SANS model
    double calcIQ(double Q) const;

    double integrand(double Q) const;

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    double ncrystalSANSSphere(double neutron_ekin) const;

    double calcXSwithIofQHelper( double neutron_ekin) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    double m_radius;
    double m_sld;
    double m_sld_solvent;
    double m_result;
    NC::IofQHelper m_Iq;//iffy! must be last for now
  };

}
#endif
