
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

#include "NCPluginFactory.hh"
#include "NCCrystallineTexture.hh"

namespace NCPluginNamespace {

  using PhysicsModel = CrystallineTexture;

  class PluginScatter final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override
    {
      return NCPLUGIN_NAME_CSTR "Model";
    }

    PluginScatter( PhysicsModel && pm ) : m_pm(std::move(pm)) {}


    NC::CrossSect crossSection(NC::CachePtr&,
                               NC::NeutronEnergy neutron_ekin,
                               const NC::NeutronDirection& ndirlab ) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection(neutron_ekin, ndirlab) };
    }


    NC::ScatterOutcome sampleScatter(NC::CachePtr&,
                                     NC::RNG& rng,
                                     NC::NeutronEnergy neutron_ekin,
                                     const NC::NeutronDirection& ndirlab) const override
    {
      return m_pm.sampleScatteringEvent( rng, neutron_ekin, ndirlab );
    }

  private:
    PhysicsModel m_pm;
  };

}

const char * NCP::PluginFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Here follows the factory logic, for how the physics model provided by the  //
// plugin should be combined with existing models in NCrystal.                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

NC::Priority
NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Must return value >0 if we should do something, and a value higher than
  //100 means that we take precedence over the standard NCrystal factory:

  if (!cfg.get_coh_elas())
    return NC::Priority::Unable;//coherent-elastic disabled, do nothing.

  //Ok, we might be applicable. Load input data and check if is something we
  //want to handle:
  if ( ! PhysicsModel::isApplicable( cfg.info() ) )
    return NC::Priority::Unable;

  //Ok, all good. Tell the framework that we want to deal with this, with a
  //higher priority than the standard factory gives (which is 100):
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr
NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:

  auto sc_pp = createStdPlaneProvider( cfg.infoPtr() );
  auto sco = cfg.createSCOrientation();
  auto sc_ourmodel
    = NC::makeSO<PluginScatter>( CrystallineTexture::createFromInfo( sco,
                                                                     cfg.info(),
                                                                     sc_pp.get() ) );

  //Now we just need to combine this with all the other physics. So ask the
  //framework to set this up, except for the component we are replacing in the
  //plugin:

  auto sc_std = globalCreateScatter( cfg.modified("coh_elas=0") );

  //Combine and return:
  return combineProcs( sc_std, sc_ourmodel );
}
