
#include "NCPluginFactory.hh"
#include "NCSansIsotropic.hh"
#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu
#include "NCSansIQCurve.hh"

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ScatterIsotropic {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    PluginScatter( SansIsotropic && pm )
      : ScatterIsotropic(NCPLUGIN_NAME_CSTR "Model"),
        m_pm(std::move(pm))
    {
    }

    double crossSectionNonOriented( double ekin ) const final
    {
      //Trivial.
      return m_pm.calcCrossSection(ekin);
    }

    void generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const final
    {
      //Specialised interface, providing delta-E and angle.
      auto outcome = m_pm.sampleScatteringEvent( *getRNG(), ekin );
      delta_ekin = outcome.ekin_final - ekin;
      angle = std::acos(outcome.mu);
    }

    void generateScattering( double ekin, const double (&indir)[3],
                             double (&outdir)[3], double& delta_ekin ) const final
    {
      //Generic interface called by most MC simulation codes.
      auto outcome = m_pm.sampleScatteringEvent( *getRNG(), ekin );
      delta_ekin = outcome.ekin_final - ekin;
      randDirectionGivenScatterMu( getRNG(), outcome.mu, indir, outdir );
    }
  protected:
    virtual ~PluginScatter() = default;
  private:
    SansIsotropic m_pm;
  };

}

const char * NCP::PluginFactory::getName() const
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

int NCP::PluginFactory::canCreateScatter( const NC::MatCfg& cfg ) const
{
  if ( ! SansIQCurve::isApplicable(*globalCreateInfo(cfg)) )
    return 0;
  return 999;
}

NC::RCHolder<const NC::Scatter> NCP::PluginFactory::createScatter( const NC::MatCfg& cfg ) const
{
  auto sc_ourmodel = NC::makeRC<PluginScatter>(SansIQCurve::createFromInfo(*globalCreateInfo(cfg), cfg.get_packfact() ));
  return sc_ourmodel;

  //fixme: to be enabled again
  // auto cfg2 = cfg.clone();
  // auto sc_std = globalCreateScatter(cfg2);
  // // Combine and return:
  // return combineScatterObjects( sc_std,
  //                               sc_ourmodel.dyncast<const NC::Scatter>() );
}
