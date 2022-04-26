
#include "NCPluginFactory.hh"
#include "NCSansIsotropic.hh"
#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu
#include "NCSansModelPicker.hh"

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

NC::Priority NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  if ( ! ( cfg.get_sans() && SansModelPicker::isApplicable(cfg.info()) ) )
    return NC::Priority::Unable;
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  auto sc_ourmodel = NC::makeSO<PluginScatter>( SansModelPicker::createFromInfo(cfg.info()) );
  auto sc_std = globalCreateScatter( cfg );
  //NOTE from TK: To get ONLY the sans component when plotting etc., you can simply put "bla.ncmat;comp=sans".
  return combineProcs( sc_std, sc_ourmodel );
}
