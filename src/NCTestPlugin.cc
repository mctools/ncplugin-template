
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCTestPlugin.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  //Note, emit all messages here and elsewhere in plugin code with NCPLUGIN_MSG
  //(or NCPLUGIN_WARN for warnings), never raw usage of std::cout or printf!

  NCPLUGIN_MSG("Testing plugin");

  //Create some test NCMAT data. For simplicity we will here base it on an
  //existing file, but add our @CUSTOMPLUGIN section (see NCPhysicsModel.cc for
  //a description of the format).:

  const std::string base_data = "stdlib::Be_sg194.ncmat";
  std::string testdata;
  {
    testdata = NC::FactImpl::createTextData(base_data)->rawDataCopy();
    std::ostringstream customsection;
    customsection << "@CUSTOM_" << pluginNameUpperCase() << "\n"
                  << "BC_pure  0  1000  100000  Gauss\n";
    testdata += customsection.str();
  }

  if ( true ) {
    NCPLUGIN_MSG("Test NCMAT data begin:");
    NCRYSTAL_RAWOUT(testdata);
    NCPLUGIN_MSG("Test NCMAT data end.");
  }

  //Let us load and exercise this testdata. We do NOT want to register the
  //testdata as a virtual file, since we want this test to not leave anything
  //registered after it is done running. So we simply put the raw data directly
  //into a MatCfg object (which is the C++ object representing a cfg-string):

  auto cfg = NC::MatCfg::createFromRawData( std::string(testdata) );
  auto cfg_ref = NC::MatCfg( base_data );
  auto scat = NC::createScatter( cfg );
  auto scat_ref = NC::createScatter( cfg_ref );

  auto wl = NC::NeutronWavelength{ 3.45 };
  auto xs_ref = scat_ref.crossSectionIsotropic( wl );
  auto xs = scat.crossSectionIsotropic( wl );

  NCPLUGIN_MSG("Ref   Be      XS @ "<<wl<<": "<<xs_ref
               <<" (should be around 9.55 barn)");
  NCPLUGIN_MSG("Dummy Be+Extn XS @ "<<wl<<": "<<xs
               <<" (should be around 5.80 barn)");

  nc_assert_always( xs_ref.dbl() > 9.4 && xs_ref.dbl() < 9.8 );
  nc_assert_always( xs.dbl() > 5.6 && xs.dbl() < 6.0 );

  //Note: We do not test the scatter sampling here, but should probably do so
  //(to test the energy changes are as expected).

  NCPLUGIN_MSG("All tests of plugin were successful!");
}
