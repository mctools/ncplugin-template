#include "NCPhysicsModel.hh"
#include <iostream>
#include "NCrystal/internal/NCMath.hh"//for NC::linspace

int main()
{
  //Very simple test which instantiates our model and calculates a few cross
  //sections and samples a few scattering events:
  //NC::VectD param{132.869, -1.33605, 0.0519763, -3.97314, 0.0510821, 5.551};
  auto pm = NCP::PhysicsModel(1, 132.869, -1.33605, 0.0519763, -3.97314, 0.0510821);
  //auto pm = NCP::PhysicsModel("teshi_xs");

  for ( auto en :  NC::logspace(-4,1,50) ) {
    std::cout << "cross section @ " << en << " eV is "
              << pm.calcCrossSection( en ) <<" barn" << std::endl;
  }

  auto rng = NC::getRNG();

  for ( auto wl : NC::linspace(9, 11, 2) ) {
    for (unsigned i = 0; i < 10; ++i) {
      auto outcome = pm.sampleScatteringEvent( *rng, NC::wl2ekin(wl) );
      if (std::isnan(outcome.mu)) {
        std::cout << "scattering @ " << wl << " Aa gives neutron with wl = "
                << NC::ekin2wl(outcome.ekin_final) <<" Aa and a scattering   "
                << std::acos(outcome.mu)*NC::kToDeg<<" degrees. The " 
                << std::endl;
      }
    }
  }



  return 0;
}
