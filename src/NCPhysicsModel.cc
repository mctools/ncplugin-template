#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCRandUtils.hh"

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_SANSND section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, it 
  //raises BadInput exceptions
  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom as
  // is usual in NCrystal):
  //
  // @CUSTOM_SANSND
  //    v                # plugin version
  //    A_1 b_1 A_2 b_2  # piecewise power law parameters
  //    Q_0 sigma_0      # boundary parameter and absolute cross section 
  //              

  //Verify we have exactly three lines and 1-4-2 numbers:
  if ( data.size() != 3 || data.at(0).size()!=1 || data.at(1).size()!=4 || data.at(2).size()!=2 )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should be three lines with 1-4-2 numbers");

  //Parse and validate values:
  double supp_version = 1.0;
  double version, A1, A2, b1, b2, Q0, sigma0 ;
  if ( ! NC::safe_str2dbl( data.at(0).at(0), version )
       || ! NC::safe_str2dbl( data.at(1).at(0), A1 )
       || ! NC::safe_str2dbl( data.at(1).at(1), b1 )
       || ! NC::safe_str2dbl( data.at(1).at(2), A2 )
       || ! NC::safe_str2dbl( data.at(1).at(3), b2 )
       || ! NC::safe_str2dbl( data.at(2).at(0), Q0 )
       || ! NC::safe_str2dbl( data.at(2).at(1), sigma0 )
       || !(Q0>0) || !(sigma0>0)) 
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" section (see the plugin readme for more info)" );
  //special warning for wrong version
  if ( ! (version==supp_version) ) 
    NCRYSTAL_THROW2( BadInput,"Invalid version specified for the "<<pluginNameUpperCase()
                     <<" plugin. Only the version "<<supp_version<<" is supported." );

  //Parsing done! Create and return our model:
  return PhysicsModel(A1, b1, A2, b2, Q0,sigma0);
}

NCP::PhysicsModel::PhysicsModel( double A1, double b1, double A2, double b2, double Q0, double sigma0 )
  : m_A1(A1),
    m_b1(b1),
    m_A2(A2),
    m_b2(b2),
    m_Q0(Q0),
    m_sigma0(sigma0)
{
  nc_assert( m_Q0 > 0.0 );
  nc_assert( m_sigma0 > 0.0 );
}

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{
  double lambda = NC::ekin2wl(neutron_ekin); //wavelength
  double k =  2*NC::kPi/lambda; //wavevector
  double total_sigma = (m_sigma0/(2*k))*(m_A1/(m_b1+2)*std::pow(m_Q0,m_b1+2) + m_A2/(m_b2+2)*std::pow(2*k,m_b2+2) - m_A2/(m_b2+2)*std::pow(m_Q0,m_b2+2));
  return total_sigma;
}

double NCP::PhysicsModel::sampleScatteringVector( NC::RandomBase& rng, double neutron_ekin ) const 
{
  double rand = rng.generate();
  double Q;
  //sample a random scattering vector Q from the inverse CDF (see plugin readme)
  double ratio_sigma = m_sigma0/calcCrossSection(neutron_ekin); //cross section over total cross section ratio
  double CDF_Q0 = (m_A1*std::pow(m_Q0, m_b1+2)/(m_b1+2))*ratio_sigma;
  if(rand < CDF_Q0){
    Q = std::pow(((m_b1+2)*rand/m_A1)/ratio_sigma, 1/m_b1+2);
  } else {
    Q = std::pow((rand/ratio_sigma - m_A1/(m_b1+2)*std::pow(m_Q0,m_b1+2) + m_A2/(m_b2+2)*std::pow(m_Q0,m_b2+2))*(m_b2+2)/m_A2, 1/m_b2+2);
  }
  return Q;
}
NCP::PhysicsModel::ScatEvent NCP::PhysicsModel::sampleScatteringEvent( NC::RandomBase& rng, double neutron_ekin ) const
{
  ScatEvent result;
  double lambda = NC::ekin2wl(neutron_ekin); //wavelength
  double k =  4*std::acos(0.0)/lambda; //wavevector

  //Implement our actual model here:
  result.ekin_final = neutron_ekin;//Elastic
  double Q = sampleScatteringVector(rng, neutron_ekin);
  result.mu = 1-0.5*(Q/k)*(Q/k);

  return result;
}


