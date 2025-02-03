#ifndef NCPlugin_CrystallineExtinction_hh
#define NCPlugin_CrystallineExtinction_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class CrystallineExtinction final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static CrystallineExtinction createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors
    
    //To account for the block size and the mosaic spread, the theories of extinction
    //introduce primary and secondary corrections to the coherent elastic neutron
    //scattering cross section for randomly oriented powders. Dependent on whether
    //blocks separated by small-angle tilts but having the same relative orientations
    //are considered correlated, two models have been developed, the correlated block
    //model and the uncorrelated block model, respectively.
    //Ref: International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616.
    //As a first investigation, the latter is studied and implemented as follows.

    CrystallineExtinction( int model_option, double l, double Gg,
                           int tilt_dist_option, double L, double c,
                           int res_option, double a0, double a1,
                           double b0, double b1,double s02,
                           double s12, double s22,
                           const NCrystal::StructureInfo& struct_info,
                           const NCrystal::HKLList& hklList );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    int m_model_option;
    double m_l;
    double m_Gg;
    int m_tilt_dist_option;
    double m_L;
	  double m_c;
    int m_res_option;
    double m_a0;
    double m_a1;
    double m_b0;
    double m_b1;
    double m_s02;
    double m_s12;
    double m_s22;
    double m_Nc;
    double m_xsectfact;
    struct DFM {
      double dspacing;
      double fsquared;
      double multiplicity;
    };
    std::vector<DFM> hklDFM;
  };

}
#endif
