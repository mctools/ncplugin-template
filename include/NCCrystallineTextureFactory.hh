#ifndef NCPlugin_CrystallineTextureFactory_hh
#define NCPlugin_CrystallineTextureFactory_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //Factory which implements logic of how the physics model provided by the
  //plugin should be combined with existing models in NCrystal:

  class CrystallineTextureFactory final : public NC::FactImpl::ScatterFactory {
  public:
    //ScatFact = Factory<FactoryType::Scatter,ProcImpl::Process,ScatterRequest>
    const char * name() const noexcept override;

    NC::Priority query( const NC::FactImpl::ScatterRequest& ) const override;
    NC::ProcImpl::ProcPtr produce( const NC::FactImpl::ScatterRequest& ) const override;
  };
}

#endif
