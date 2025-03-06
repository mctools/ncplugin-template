
//Define exported symbols used by NCrystal to activate the plugin. This is
//mostly generic, except for the registerPlugin function which can potentially
//be modified in some use-cases.

#define NCPLUGIN_BOILERPLATE_CC
#include "NCrystal/NCPluginBoilerplate.hh"

#include "NCPluginFactory.hh"
#include "NCTestPlugin.hh"

void NCP::registerPlugin()
{
  //This function is required for the plugin to work. It should register
  //factories, and potentially other stuff as appropriate for the plugin (like
  //adding in-mem data files, adding test functions, ...).
  NC::FactImpl::registerFactory(std::make_unique<NCP::PluginFactory>());
  NC::Plugins::registerPluginTestFunction( std::string("test_") + pluginName(),
                                           customPluginTest );
};
