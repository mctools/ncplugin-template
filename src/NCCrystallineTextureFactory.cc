
#include "NCCrystallineTextureFactory.hh"
#include "NCCrystallineTexture.hh"
//#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "Model"; }
    PluginScatter( CrystallineTexture && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect crossSection(NC::CachePtr&, NC::NeutronEnergy neutron_ekin, const NC::NeutronDirection& ndirlab ) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection(neutron_ekin, ndirlab) };
    }

    NC::ScatterOutcome sampleScatter(NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy neutron_ekin, const NC::NeutronDirection& ndirlab) const override
    {
      return m_pm.sampleScatteringEvent( rng, neutron_ekin, ndirlab );
    }

  private:
    CrystallineTexture m_pm;
  };

}

const char * NCP::CrystallineTextureFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// Here follows the factory logic, for how the physics model provided by the    //
// plugin should be combined with existing models in NCrystal.                  //
//                                                                              //
// In the silly example here, we want our custom physics model to replace the   //
// existing incoherent-elastic model of NCrystal with our own model.            //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

NC::Priority NCP::CrystallineTextureFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Must return value >0 if we should do something, and a value higher than
  //100 means that we take precedence over the standard NCrystal factory:
  //if (!cfg.get_coh_elas())
  //  return NC::Priority::Unable;//coherent-elastic disabled, never do anything.

  //Ok, we might be applicable. Load input data and check if is something we
  //want to handle:

  if ( ! CrystallineTexture::isApplicable( cfg.info() ) )
    return NC::Priority::Unable;

  //Ok, all good. Tell the framework that we want to deal with this, with a
  //higher priority than the standard factory gives (which is 100):
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr NCP::CrystallineTextureFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:

  auto sc_pp = createStdPlaneProvider( cfg.infoPtr() );
  NC::SCOrientation sco = cfg.createSCOrientation();
  auto sc_ourmodel = NC::makeSO<PluginScatter>( CrystallineTexture::createFromInfo( sco, cfg.info(), sc_pp.get() ) );

  //Now we just need to combine this with all the other physics
  //(i.e. Bragg+inelastic).  So ask the framework to set this up, except for
  //incoherent-elastic physics of course since we are now dealing with that
  //ourselves:

  auto sc_std = globalCreateScatter( cfg.modified("coh_elas=0") );

  //Combine and return:
  return combineProcs( sc_std, sc_ourmodel );
}
