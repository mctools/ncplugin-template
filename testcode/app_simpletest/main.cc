#include "NCPhysicsModel.hh"
#include <iostream>
#include "NCrystal/internal/NCMath.hh"//for NC::linspace
#include "NCExtraTestUtils.hh"

int main()
{
  //Very simple test which instantiates our model and calculates a few cross
  //sections and samples a few scattering events:
  //NC::VectD param{132.869, -1.33605, 0.0519763, -3.97314, 0.0510821, 5.551};
  //auto pm = NCP::PhysicsModel(NCrystalPlugin_SANSND::PhysicsModel::Model::PPF, 132.869, 1.33605, 0.0519763, 3.97314, 0.0510821);
  //NC::VectD param{50};
  NC::VectD param{132.869, 1.33605, 0.0519763, 3.97314, 0.0510821, 1};
  auto pm = NCP::PhysicsModel(NCrystalPlugin_SANSND::PhysicsModel::Model::PPF, param);

  /*for ( auto en :  NC::logspace(-4,1,50) ) {
    std::cout << "cross section @ " << en << " eV is "
              << pm.calcCrossSection( en ) <<" barn" << std::endl;
  }*/

  //auto rng = NC::getRNG();

  /*for ( auto wl : NC::linspace(9, 11, 2) ) {
    for (unsigned i = 0; i < 10; ++i) {
      auto outcome = pm.sampleScatteringEvent( *rng, NC::wl2ekin(wl) );
      if (std::isnan(outcome.mu)) {
        std::cout << "scattering @ " << wl << " Aa gives neutron with wl = "
                << NC::ekin2wl(outcome.ekin_final) <<" Aa and a scattering   "
                << std::acos(outcome.mu)*NC::kToDeg<<" degrees. The " 
                << std::endl;
      }
    }
  }*/
  double avg;
  auto angles = NCPluginTestCode::sampleAngles(pm, NC::wl2ekin(10), 1000);
  for (size_t i = 0; i < 1000; i++)
  {
    avg += angles[i]*180/(2*3.1421)/1000;
  }
    std::cout << avg << std::endl;

  auto ekins = NC::logspace(-5,1,1000);
  auto xs = NCPluginTestCode::sampleXS(pm, ekins, 1000);
  std::cout << xs[0] << std::endl;
  //for (size_t i = 0; i < ekins.size(); i++)
  //{
  //  std::cout << xs[i] << std::endl;
  //} 
  

  return 0;
}
