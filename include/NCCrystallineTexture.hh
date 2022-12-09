#ifndef NCPlugin_CrystallineTexture_hh
#define NCPlugin_CrystallineTexture_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCMatrix.hh"

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
    static CrystallineTexture createFromInfo( const NC::Info& );//will raise BadInput in case of syntax errors

    //The crystalline texuture or preferred orientation correction is taken into account by introducing 
    //the cylindrically symmetric Pole-Density Distribution Function (PDDF) P_hkl(theta_hkl) which 
    //depends on the orientation angle theta_hkl, i.e., angle between the preferred orientation and
    //the plan vectors tau_hkl.
    //Various types of PDDFs are reported in the literature.
    //As a first try, we investigate the March-Dollase model.

    //Constructor gets constant cross section value, and the neutron wavelength
    //cutoff:
    CrystallineTexture( NCrystal::Vector& preferred_orientation1, double p1, double f1,
                        NCrystal::Vector& preferred_orientation2, double p2, double f2,
                        const NCrystal::StructureInfo& struct_info,
                        const NCrystal::HKLList& hkl_list );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    NCrystal::Vector& m_preferred_orientation1;
    double m_p1;
    double m_f1;
    NCrystal::Vector& m_preferred_orientation2;
    double m_p2;
    double m_f2;
    const NCrystal::StructureInfo& m_struct_info;
    const NCrystal::HKLList& m_hkl_list;
  };

}
#endif
