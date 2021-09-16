
#include "NCPluginFactory.hh"
#include "NCSansIsotropic.hh"
#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu
#include "NCSansIQCurve.hh"

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "Model"; }
    PluginScatter( SansIsotropic && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection(ekin.dbl()) };
    }

    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const override
    {
      auto outcome = m_pm.sampleScatteringEvent( rng, ekin.dbl() );
      return { NC::NeutronEnergy{outcome.ekin_final}, NC::CosineScatAngle{outcome.mu} };
    }

  private:
    SansIsotropic m_pm;
  };

}

const char * NCP::PluginFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

NC::Priority NCP::PluginFactory::query( const NC::MatCfg& cfg ) const
{
  if ( ! SansIQCurve::isApplicable(*globalCreateInfo(cfg)) )
    return NC::Priority::Unable;
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr NCP::PluginFactory::produce( const NC::MatCfg& cfg ) const
{
  // cfg.get_packfact()
  auto sc_ourmodel = NC::makeSO<PluginScatter>(SansIQCurve::createFromInfo(globalCreateInfo(cfg) ));
  return sc_ourmodel;

  //fixme: to be enabled again
  // auto cfg2 = cfg.clone();
  // auto sc_std = globalCreateScatter(cfg2);
  //[Comment from TK: No need to clone cfg as cfg2 if you are not going to change any parameters.]
  //
  // // Combine and return:
  // return combineProcs( sc_std, sc_ourmodel );
}
