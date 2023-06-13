
#include "NCPluginFactory.hh"
#include "NCrystal/internal/NCRandUtils.hh" // for randDirectionGivenScatterMu
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCSANSDiluteIsotropicScatter.hh"
#include "NCrystal/NCFactImpl.hh"

namespace NCPluginNamespace {
  namespace {

    double sas_3j1x_x(double q)
    {
      //NB: q is actually q*sphere_radius
      constexpr double SPH_J1C_CUTOFF = 0.1;
      if ( NC::ncabs(q) < SPH_J1C_CUTOFF ) {
        const double q2 = q*q;
        return (1.0 + q2*(-3./30. + q2*(3./840. + q2*(-3./45360.))));// + q2*(3./3991680.)))));
      } else {
        double sin_q, cos_q;
        NC::sincos(q,cos_q, sin_q);
        return 3.0*(sin_q/q - cos_q)/(q*q);
      }
    }
}



  // class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  // public:

  //   //The factory wraps our custom PhysicsModel helper class in an NCrystal API
  //   //Scatter class.

  //   const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "Model"; }
  //   PluginScatter( PhysicsModel && pm ) : m_pm(std::move(pm)) {}

  //   NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin) const override
  //   {
  //     return NC::CrossSect{ m_pm.calcCrossSection(ekin.dbl()) };
  //   }

  //   NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const override
  //   {
  //     auto outcome = m_pm.sampleScatteringEvent( rng, ekin.dbl() );
  //     return { NC::NeutronEnergy{outcome.ekin_final}, NC::CosineScatAngle{outcome.mu} };
  //   }

  // private:
  //   PhysicsModel m_pm;
  // };

}

const char * NCP::PluginFactory::name() const noexcept
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

NC::Priority NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& request ) const
{
  return ( ( request.get_sans() && hasCustomDataForSANSPlugin( request.info(), "HARDSPHERESANS" ) )
           ? Priority::OnlyOnExplicitRequest
           : Priority::Unable );
}

NC::ProcImpl::ProcPtr NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& request ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:
  nc_assert( request.get_sans() );
  auto datas = extractCustomDataForSANSPlugin( request.info(), "HARDSPHERESANS" );
  nc_assert( !datas.empty() );

  NC::ProcImpl::ProcComposition::ComponentList process_list;

  //First add all the physics provided by other code:
  process_list.emplace_back( globalCreateScatter(request.modified("scatfactory=!hardspheresans")) );//todo: make more robust

  //Now add our SANS model as appropriate (this might be more than one
  //instance in case of complicated material setups):
  for ( auto& data : datas ) {

    //Decode and check syntax in custom section data:
    const auto& lines = data.customData;
    double sphere_radius;
    if ( lines.size()!=1
         || lines.at(0).size()!=1
         || !NC::safe_str2dbl(lines.at(0).at(0),sphere_radius)
         || !(sphere_radius>0.0)
         || !(sphere_radius<1e6) ) {
      NCRYSTAL_THROW2(BadInput,"Syntax error in @CUSTOM_"<< "HARDSPHERESANS"
                      <<" section. Expects a single positive number (the sphere radius in Angstrom).");
    }

    //Our SANS process
    SANSDiluteIsotropicModel model;
    model.pofq = [sphere_radius]( double q ) {
      return NC::ncsquare( sas_3j1x_x(q*sphere_radius) );
    };
    model.form_volume = NC::nccube(sphere_radius)*(4.0/3.0)*NC::kPi;
    process_list.emplace_back( NC::makeSO<SANSDiluteIsotropicScatter>( data.scale, model ) );
  }

  return NC::ProcImpl::ProcComposition::consumeAndCombine(std::move(process_list));
}
