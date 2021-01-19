#include "NCSansIQCurve.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include <vector>

bool NCP::SansIQCurve::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::SansIsotropic NCP::SansIQCurve::createFromInfo( const NC::Info& info, double packfact)
{
  if(packfact<0. || packfact>1.)
    NCRYSTAL_THROW2(BadInput, "Packing factor in @CUSTOM_"<<pluginNameUpperCase()
                    <<" should be in the range between (0,1], but " << packfact <<" is given.");
  auto iq = SansIQCurve(info, packfact);
  return SansIsotropic(iq.getQ(), iq.getI());
}


bool NCP::SansIQCurve::calSDL(const NC::Info& info, double &scatLenDensity, double &numberDensity) const
{
  scatLenDensity= 0.;
  numberDensity =0.;
  //calculate scattering length density from the dynamic info
  if(info.hasDynamicInfo()&&info.hasNumberDensity())
  {
    numberDensity = info.getNumberDensity()*m_densityScale;   // in atoms/Aa^3
    for (auto& dyn : info.getDynamicInfoList())
    {
      double scl = dyn->atomDataSP()->coherentScatLen(); //in sqrt(barn)
      double frac = dyn->fraction();
      scatLenDensity += scl*frac*numberDensity;
    }
  }
  else if(info.hasStructureInfo()&&info.hasAtomPositions())
  {
    auto &strInfo = info.getStructureInfo();
    double perVolume = 1./strInfo.volume*m_densityScale;//Aa^3

    for(auto it = info.atomInfoBegin(); it != info.atomInfoEnd(); ++it)
    {
      double scl = it->data().coherentScatLen(); //in sqrt(barn)
      scatLenDensity += scl*perVolume*it->number_per_unit_cell;
      numberDensity += it->number_per_unit_cell*perVolume;
    }
  }
  else
    return false;
  return true;
}

NCP::SansIQCurve::SansIQCurve( const NC::Info& info, double packfact )
:m_densityScale(packfact)
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  double sdl(0), density(0);
  calSDL(info, sdl, density);
  printf("sdl %g, density %g\n",sdl, density);

  //Verify we have two lines and have identical number of elements
  if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  if(data.at(0).at(0)!= "Q" || data.at(1).at(0)!= "I")
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                  <<" the first and second line should be I and Q, respectively.");

  //Parse and validate values:
  m_I.reserve(data.at(0).size());
  m_Q.reserve(data.at(0).size());

  for(unsigned i=1;i<data.at(0).size();i++)
  {
    double temp(0.);
    if ( !NC::safe_str2dbl( data.at(1).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section I" );
    m_I.push_back(temp);

    if ( !NC::safe_str2dbl( data.at(0).at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section Q" );
    m_Q.push_back(temp);
  }
}
