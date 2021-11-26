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
    //after debugging
    //enum class model_def : unsigned { FILE=0, PPF=1, GPF=2, HSFBA=3 };
    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static PhysicsModel createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors

    //Constructor gets the filename of the input data file:
    PhysicsModel(std::string filename, double norm = 1 );
    //Constructor gets the models string and the param:
    PhysicsModel( int model, double p0, double p1, double p2, double p3, double p4 );
    //Constructor gets the models string and nanoparticle radius R :
    PhysicsModel( int model, double mono_R);
    //Constructor gets the models string and the dist R file:
    PhysicsModel( int model, std::string filename );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;
    //Sample scattering vector from inverse CDF (rng is random number stream).
    double sampleScatteringVector( NC::RNG& rng, double neutron_ekin ) const;

    //Sample scattering event. Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    int m_model;
    NC::Optional<NC::VectD> m_param;
    NC::Optional<NC::IofQHelper> m_helper;
    //for theoretical NP_FBA
    NC::Optional<double> m_mono_R;
    //NC::Optional<std::pair<NC::VectD,NC::VectD>> m_R_freq;
    //for reading from file
    NC::Optional<double> m_norm;
  };

}
#endif
