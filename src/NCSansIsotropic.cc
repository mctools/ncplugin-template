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

NCP::SansIsotropic::SansIsotropic( const std::vector<double>& x, const std::vector<double>& f )
  : m_Iq(std::make_unique<LookUpTable>(x, f))
{
  m_Iq->sanityCheck();
}

double NCP::SansIsotropic::calcCrossSection( double neutron_ekin ) const
{
  //fixme
  return 10.0+neutron_ekin*0;
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
