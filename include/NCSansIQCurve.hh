#ifndef NCPlugin_SansIQCurve_hh
#define NCPlugin_SansIQCurve_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

#include "NCSansIsotropic.hh"
#include <string>

namespace NCPluginNamespace {

  //This class helps to create an I(Q) curve from info object for the Isotropic sans model

  class SansIQCurve final : public NC::MoveOnly {
  public:
    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static SansIsotropic createFromInfo( const NC::Info& ); //will raise BadInput in case of syntax errors

  private:
    SansIQCurve( const NC::Info& );
    bool calSDL(const NC::Info& info, double &scatLenDensity, double &numberDensity) const;
    // not have euough info to calculate material coherent scattering length density and number density
    const std::vector<double>& getQ() const {return m_Q;}
    const std::vector<double>& getI() const {return m_I;}


  private:
    //Data members:
    std::string m_solvantCfg;
    std::vector<double> m_I, m_Q;

  };
}
#endif
