#include "NCSansIsotropic.hh"
//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include <vector>

bool NCP::SansIsotropic::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::SansIsotropic NCP::SansIsotropic::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <sigmavalue> <wavelength threshold value>
  //

  //Verify we have two lines and have identical number of elements
  if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  if(data.at(0).at(0)!= "I" || data.at(1).at(0)!= "Q")
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                  <<" the first and second line should be I and Q, respectively.");

  //Parse and validate values:
  std::vector<double> I, Q;
  I.reserve(data.at(0).size());
  Q.reserve(data.at(0).size());

  for(unsigned i=1;i<data.at(0).size();i++)
  {
    double temp(0.);
    if ( !NC::safe_str2dbl( data.at(0).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section I" );
    I.push_back(temp);

    if ( !NC::safe_str2dbl( data.at(1).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section S" );
    Q.push_back(temp);
  }

  //Parsing done! Create and return our model:
  return SansIsotropic(I, Q);
}

double NCP::SansIsotropic::ekin2k(double ekin) const
{
  return 2.*NC::kPi/NC::ekin2wl(ekin);
}


NCP::SansIsotropic::SansIsotropic( const std::vector<double>& Q, const std::vector<double>& intensity )
  : m_Iq(std::make_unique<LookUpTable>(Q, intensity, NCP::LookUpTable::Extrapolate::kConst_Zero)),
   m_xs(std::make_unique<LookUpTable>())
{
  m_Iq->sanityCheck();
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
    double lastK=ekin2k(envec.back());
    accumIntegrand = xsvec.back()*lastK*lastK;
    accumIntegrand += NC::kPi*(Q[i]-Q[i-1])*(intensity[i]*Q[i]+intensity[i-1]*Q[i-1]);
    double k = Q[i]*0.5;
    envec.push_back(NC::wl2ekin(2*NC::kPi/k));
    xsvec.push_back(accumIntegrand/(k*k));
  }

  m_xs=std::make_unique<LookUpTable>(envec, xsvec, NCP::LookUpTable::Extrapolate::kZero_Zero);
}

double NCP::SansIsotropic::calcCrossSection( double neutron_ekin ) const
{
  //fixme
  return m_xs->get(neutron_ekin);
}

NCP::SansIsotropic::ScatEvent NCP::SansIsotropic::sampleScatteringEvent( NC::RandomBase& rng, double neutron_ekin ) const
{
  ScatEvent result;

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
