#ifndef NCPlugin_SansHelper_hh
#define NCPlugin_SansHelper_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include <string>

namespace NCPluginNamespace {

  bool calSDL(const NC::Info& info, double &scatLenDensity, double &numberDensity);
  NC::Info::CustomSectionData::const_iterator findCustomLineIter(const NC::Info::CustomSectionData& data, const std::string& keyword, bool check=true);

}
#endif
