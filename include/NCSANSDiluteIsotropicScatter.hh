#ifndef NCPlugin_SANSGenericFct_hh
#define NCPlugin_SANSGenericFct_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/internal/NCIofQHelper.hh"
#include "NCrystal/internal/NCSANSUtils.hh"//For SANSScaleFactor

namespace NCPluginNamespace {

  struct SANSDiluteIsotropicModel {
    using PofQFct = std::function<double(double)>;
    PofQFct pofq;
    double form_volume = 0.0;
    //TODO: std::string name;
    //TODO: qgrid...
  };

  class SANSDiluteIsotropicScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "SANSDiluteIsotropicScatter"; }
    SANSDiluteIsotropicScatter( NC::SANSScaleFactor, SANSDiluteIsotropicModel );

    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin) const override;
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic( NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const override;

  private:
    double m_xsscale;
    NC::IofQHelper m_iofq;
  };

}

#endif
