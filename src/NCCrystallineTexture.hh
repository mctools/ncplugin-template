#ifndef NCPlugin_CrystallineTexture_hh
#define NCPlugin_CrystallineTexture_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"

namespace NCPluginNamespace {

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class CrystallineTexture final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static CrystallineTexture createFromInfo( const NC::SCOrientation&, const NC::Info&, NC::PlaneProvider * = nullptr );//will raise BadInput in case of syntax errors

    //The crystalline texuture or preferred orientation correction is taken into account by introducing 
    //the cylindrically symmetric Pole-Density Distribution Function (PDDF) or preferred orientation
    //distribution function P_hkl(lambda, theta_hkl) which depends on the orientation angle theta_hkl,
    //i.e., angle between the preferred orientation and the plan vectors hkl.
    //Various types of PDDFs are reported in the literature.
    //As a first try, we investigate the modified March-Dollase model (Sato 2011).

    //Constructor gets constant cross section value, and the neutron wavelength
    //cutoff:
    CrystallineTexture( const NC::SCOrientation&,
                        const NCrystal::Vector& preferred_orientation1, double R1, double f1,
                        const NCrystal::Vector& preferred_orientation2, double R2, double f2,
                        const NCrystal::StructureInfo& struct_info,
                        NC::PlaneProvider * std_pp = nullptr );

    //Provide cross sections for a given neutron:
    double calcCrossSection( NC::NeutronEnergy, const NC::NeutronDirection&) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).

    NC::ScatterOutcome sampleScatteringEvent( NC::RNG&, NC::NeutronEnergy, const NC::NeutronDirection& ) const;

  private:
    //Data members:
    NCrystal::Vector m_preferred_orientation1;
    double m_R1;
    double m_f1;
    NCrystal::Vector m_preferred_orientation2;
    double m_R2;
    double m_f2;
    struct HKLPlane {
      NCrystal::Vector hkl; //vector (h,k,l)
      double d_hkl; //dspacing
      double strength; //dspacing*fsq*xsectfact
    };
    std::vector<HKLPlane> m_hklPlanes;
    NCrystal::RotMatrix m_lab2cry;
    NCrystal::RotMatrix m_reclat;

  };

}
#endif
