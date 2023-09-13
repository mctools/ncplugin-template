#include "NCExtraTestUtils.hh"

std::vector<double> NCPluginTestCode::sampleAngles( const NCP::PhysicsModel& pm, double ekin, unsigned nvalues )
{
  std::vector<double> angles;
  angles.reserve(nvalues);
  auto rng = NC::getRNG();
  for (unsigned i = 0; i < nvalues; ++i)
    angles.push_back( std::acos( pm.sampleScatteringEvent(rng,ekin).mu ) );
  return angles;
}

std::vector<double> NCPluginTestCode::sampleXS( const NCP::PhysicsModel& pm, std::vector<double> ekins, unsigned nvalues )
{
  std::vector<double> xs;
  xs.reserve(nvalues);
  auto rng = NC::getRNG();
  for (unsigned i = 0; i < nvalues; ++i)
    xs.push_back( pm.calcCrossSection(ekins.at(i)) );
  return xs;
}
