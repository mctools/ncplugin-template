#include "NCPhysicsModel.hh"
#include <iostream>
#include "NCrystal/internal/NCMath.hh"//for NC::linspace

int main()
{
  //Very simple test which instantiates our model and calculates a few cross
  //sections and samples a few scattering events:

  auto pm = NCP::PhysicsModel( 2.0/*sigma below lambda_cutoff*/, 5.0/*lambda cutoff*/ );

  for ( auto wl : NC::linspace(0.01, 8.0, 20) ) {
    std::cout << "cross section @ " << wl << " Aa is "
              << pm.calcCrossSection( NC::wl2ekin(wl) ) <<" barn" << std::endl;
  }

  auto rng = NC::getRNG();

  for ( auto wl : NC::linspace(1.0, 8.0, 3) ) {
    for (unsigned i = 0; i < 10; ++i) {
      auto outcome = pm.sampleScatteringEvent( *rng, NC::wl2ekin(wl) );
      std::cout << "scattering @ " << wl << " Aa gives neutron with wl = "
                << NC::ekin2wl(outcome.ekin_final) <<" Aa and a scattering angle of "
                << std::acos(outcome.mu)*NC::kToDeg<<" degrees" << std::endl;

    }
  }



  return 0;
}
