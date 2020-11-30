
#include "NCPluginFactory.hh"
#include "NCPhysicsModel.hh"
#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ScatterIsotropic {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    PluginScatter( PhysicsModel && pm )
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
    PhysicsModel m_pm;
  };

}

const char * NCP::PluginFactory::getName() const
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

int NCP::PluginFactory::canCreateScatter( const NC::MatCfg& cfg ) const
{
  //Must return value >0 if we should do something, and a value higher than
  //100 means that we take precedence over the standard NCrystal factory:
  if (!cfg.get_incoh_elas())
    return 0;//incoherent-elastic disabled, never do anything.

  //Ok, we might be applicable. Load input data and check if is something we
  //want to handle:

  if ( ! PhysicsModel::isApplicable(*globalCreateInfo(cfg)) )
    return 0;

  //Ok, all good. Tell the framework that we want to deal with this, with a
  //higher priority than the standard factory gives (which is 100):
  return 999;
}

NC::RCHolder<const NC::Scatter> NCP::PluginFactory::createScatter( const NC::MatCfg& cfg ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:

  auto sc_ourmodel = NC::makeRC<PluginScatter>(PhysicsModel::createFromInfo(*globalCreateInfo(cfg)));

  //Now we just need to combine this with all the other physics
  //(i.e. Bragg+inelastic).  So ask the framework to set this up, except for
  //incoherent-elastic physics of course since we are now dealing with that
  //ourselves:
  auto cfg2 = cfg.clone();
  cfg2.set_incoh_elas(false);
  auto sc_std = globalCreateScatter(cfg2);

  //Combine and return:
  return combineScatterObjects( sc_std,
                                sc_ourmodel.dyncast<const NC::Scatter>() );
}
