#include "NCExtraTestUtils.hh"

std::vector<double> NCPluginTestCode::sampleAngles( const NCP::SansIsotropic& pm, double ekin, unsigned nvalues )
{
  std::vector<double> angles;
  angles.reserve(nvalues);
  auto rng = NC::defaultRNG();
  for (unsigned i = 0; i < nvalues; ++i)
    angles.push_back( std::acos( pm.sampleScatteringEvent(*rng,ekin).mu ) );
  return angles;
}
