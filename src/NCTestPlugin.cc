
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
#include "NCrystal/dump/NCDump.hh"

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  NCRYSTAL_MSG("Testing plugin "<<pluginName());

  //Create some test NCMAT data. For simplicity we will here base it on an
  //existing file, but add our @CUSTOMPLUGIN section (see NCPhysicsModel.cc for
  //a description of the format).:

  //Make sure we can load:
  NCRYSTAL_MSG("Loading 'plugins::SANSND/ncplugin-SANSND_nanodiamond.ncmat'");
  NCRYSTAL_MSG("  -> info object:");
  auto info
    = NC::createInfo("plugins::SANSND/ncplugin-SANSND_nanodiamond.ncmat");
  NC::dump(info);

  NCRYSTAL_MSG("  -> scatter process:");
  auto sc
    = NC::createScatter("plugins::SANSND/ncplugin-SANSND_nanodiamond.ncmat");
  NCRYSTAL_RAWOUT(sc.underlying().jsonDescription())

  //Sanity check, some standard material does NOT have sans:
  nc_assert_always( NC::createScatter("stdlib::Al_sg225.ncmat"
                                      ";comp=sans").isNull() );
  //But our nanodiamonds do:
  nc_assert_always( ! NC::createScatter("plugins::SANSND/ncplugin-SANSND_nanodiamond.ncmat"
                                        ";comp=sans").isNull() );


  NCRYSTAL_MSG("All tests of plugin "<<pluginName()<<" were successful!");
}
