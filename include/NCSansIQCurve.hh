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
    enum IqCalType {
      kDirectLoad,
      kHardSphere,
      kUndefined,
    };
  public:
    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):
    static bool isApplicable( const NC::Info& );
    static SansIsotropic createFromInfo( NCrystal::shared_obj<const NCrystal::Info> info ); //will raise BadInput in case of syntax errors
    bool calSDL(const NC::Info& info, double &scatLenDensity, double &numberDensity) const;

  private:
    SansIQCurve(const NC::Info& info);
    IqCalType getIqCalType(const NC::Info::CustomSectionData& data) const;
    void IqDirectLoad(const NC::Info::CustomSectionData& data);
    void IqHardSphere(const NC::Info::CustomSectionData& data, const NC::Info& info);
    const std::vector<double>& getQ() const {return m_Q;}
    const std::vector<double>& getI() const {return m_I;}

  private:
    //Data members:
    double m_packfact, m_volfact; //fixme: this two parameters are set to unity for the time being
    std::string m_solvantCfg;
    std::vector<double> m_I, m_Q;
  };
}
#endif
