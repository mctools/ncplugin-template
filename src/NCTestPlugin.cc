
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  NCRYSTAL_MSG("Testing plugin "<<pluginName());

  nc_assert_always( pluginName() == "CrysText" );//sanity check

  nc_assert_always( NC::FactImpl::hasScatterFactory("CrysTextFactory" ) );
  NCRYSTAL_MSG("Verified presence of scatter factory named 'CrysTextFactory'");

  {
    NCRYSTAL_MSG("Loading 'plugins::CrysText/Al_sg225.ncmat'");
    auto info = NC::createInfo("plugins::CrysText/Al_sg225.ncmat");
    NCRYSTAL_MSG("  -> info object:");
    NC::dump(info);
    auto sc = NC::createScatter("plugins::CrysText/Al_sg225.ncmat"
                                ";dir1=@crys_hkl:0,1,0@lab:0,1,0"
                                ";dir2=@crys_hkl:1,0,0@lab:1,0,0"
                                ";mos=0.1deg"
                                );
    NCRYSTAL_MSG("  -> scatter process:");
    NCRYSTAL_RAWOUT(sc.underlying().jsonDescription())
  }

  //TODO: More comprehensive sanity checks here!!!!

  NCRYSTAL_MSG("All tests of plugin "<<pluginName()<<" were successful!");
}
