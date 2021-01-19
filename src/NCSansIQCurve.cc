#include "NCSansIQCurve.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include <vector>

bool NCP::SansIQCurve::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::SansIsotropic NCP::SansIQCurve::createFromInfo( const NC::Info& info)
{
  auto iq = SansIQCurve(info);
  return SansIsotropic(iq.getQ(), iq.getI());
}


NCP::SansIQCurve::SansIQCurve( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );


  //Verify we have two lines and have identical number of elements
  if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  if(data.at(0).at(0)!= "Q" || data.at(1).at(0)!= "I")
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                  <<" the first and second line should be I and Q, respectively.");

  //Parse and validate values:
  m_I.reserve(data.at(0).size());
  m_I.reserve(data.at(0).size());

  for(unsigned i=1;i<data.at(0).size();i++)
  {
    double temp(0.);
    if ( !NC::safe_str2dbl( data.at(1).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section I" );
    m_I.push_back(temp);

    if ( !NC::safe_str2dbl( data.at(0).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section S" );
    m_I.push_back(temp);
  }
}
