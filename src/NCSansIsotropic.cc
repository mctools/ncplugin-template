#include "NCSansIsotropic.hh"
//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include <vector>


NCP::SansIsotropic::SansIsotropic( const std::vector<double>& Q, const std::vector<double>& intensity )
: m_xs(std::unique_ptr<LookUpTable>(nullptr)),
  m_sansQDist(std::unique_ptr<NC::PointwiseDist>(nullptr))
{
  // differential cross section
  std::vector<double> diffXS;
  diffXS.reserve(Q.size());
  for(unsigned i=0;i<Q.size();i++)
  {
    diffXS.push_back(Q[i]*intensity[i]);
  }
  m_sansQDist = std::make_unique<NC::PointwiseDist>(Q, diffXS);

  // total xs
  std::vector<double> envec, xsvec;
  envec.reserve(Q.size());
  xsvec.reserve(Q.size());

  //When Q approching zero, I(Q) is a constant,
  //the cross section in this case is 4*pi*I(Q).
  //This is how the first cross section point is calculated

  double en=NC::wl2ekin(4*NC::kPi/Q[0]);
  envec.push_back(en);
  xsvec.push_back(4*NC::kPi*intensity[0]);

  double accumIntegrand(0.);

  for(size_t i=1;i<Q.size();i++)
  {
    double lastK=2.*NC::kPi/NC::ekin2wl(envec.back()); //ekin to ki
    accumIntegrand = xsvec.back()*lastK*lastK;
    accumIntegrand += NC::kPi*(Q[i]-Q[i-1])*(intensity[i]*Q[i]+intensity[i-1]*Q[i-1]);
    double k = Q[i]*0.5;
    envec.push_back(NC::wl2ekin(2*NC::kPi/k));
    xsvec.push_back(accumIntegrand/(k*k));
  }
  m_xs=std::make_unique<LookUpTable>(envec, xsvec, NCP::LookUpTable::Extrapolate::kConst_OverSqrtX);
}

double NCP::SansIsotropic::calcCrossSection( double neutron_ekin ) const
{
  if(m_xs.get())
    return m_xs->get(neutron_ekin);
  else
    return 0.;
}

NCP::SansIsotropic::ScatEvent NCP::SansIsotropic::sampleScatteringEvent( NC::RandomBase& rng, double neutron_ekin ) const
{
  ScatEvent result;
  result.ekin_final = neutron_ekin;

  double kappa = 2*NC::kPi/NC::ekin2wl(neutron_ekin);
  double Q(1e10);
  while (Q>2*kappa) {
    Q = m_sansQDist->percentile(rng.generate());
  }
  result.mu = 1-Q*Q/(2*kappa*kappa);

  //fixme

  // if ( ! (neutron_ekin > m_cutoffekin) ) {
  //   //Special case: We are asked to sample a scattering event for a neutron
  //   //energy where we have zero cross section! Although in a real simulation we
  //   //would usually not expect this to happen, users with custom code might
  //   //still generate such calls. The only consistent thing to do when the cross
  //   //section is zero is to not change the neutron state parameters, which means:
  //   result.ekin_final = neutron_ekin;
  //   result.mu = 1.0;
  //   return result;
  // }
  //
  // //Implement our actual model here. Of course it is trivial for the example
  // //model. For a more realistic or complicated model, it might be that
  // //additional helper classes or functions should be created and used, in order
  // //to keep the code here manageable:
  //
  // result.ekin_final = neutron_ekin;//Elastic
  // result.mu = randIsotropicScatterMu(&rng);//Isotropic.

  return result;
}
