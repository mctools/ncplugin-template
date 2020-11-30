#ifndef NCPlugin_Factory_hh
#define NCPlugin_Factory_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //Factory which implements logic of how the physics model provided by the
  //plugin should be combined with existing models in NCrystal:

  class PluginFactory final : public NC::FactoryBase {
  public:
    const char * getName() const final;
    int canCreateScatter( const NC::MatCfg& ) const final;
    NC::RCHolder<const NC::Scatter> createScatter( const NC::MatCfg& ) const final;
  };
}

#endif
