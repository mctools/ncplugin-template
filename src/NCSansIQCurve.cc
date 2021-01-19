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

NCP::SansIQCurve::IqCalType NCP::SansIQCurve::getIqCalType(const NC::Info::CustomSectionData& data) const
{
  for(auto line:data)
  {
    if(line.at(0)=="DirectLoad")
      return IqCalType::kDirectLoad;
  }
  return IqCalType::kUndefined;
}


void NCP::SansIQCurve::IqDirectLoad(const NC::Info::CustomSectionData& data)
{
  //Verify we have 3 lines and 2 vectors has identical number of elements
  // if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
  //   NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
  //                   <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  NC::Info::CustomSectionData::const_iterator dataQ(data.end()), dataI(data.end());
  for(auto line=data.begin();line!=data.end();++line)
  {
    if(line->at(0)=="Q")
      dataQ=line;
    else if(line->at(0)=="I")
      dataI=line;
  }

  if(dataQ==data.end() || dataI==data.end())
    NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" direct load, do not found both I and Q vector");

  //Parse and validate values:
  m_I.reserve(dataQ->size());
  m_Q.reserve(dataQ->size());

  for(unsigned i=1;i<dataQ->size();i++)
  {
    double temp(0.);
    if ( !NC::safe_str2dbl( dataI->at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section I" );
    m_I.push_back(temp);

    if ( !NC::safe_str2dbl( dataQ->at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section Q" );
    m_Q.push_back(temp);
  }
}

NCP::SansIQCurve::SansIQCurve( const NC::Info& info, double packfact )
:m_densityScale(packfact)
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.


  unsigned numSec =  info.countCustomSections( pluginNameUpperCase() );

  for(unsigned si=0;si<numSec;si++)
  {
    NC::Info::CustomSectionData data = info.getCustomSection( pluginNameUpperCase(), si );

    switch(getIqCalType(data)) {
      case kDirectLoad:
        IqDirectLoad(data);
        break;
      default :
        NCRYSTAL_THROW2(BadInput," @CUSTOM_"<<pluginNameUpperCase()<< " with undefined load method");
    }

    double sdl(0), density(0);
    calSDL(info, sdl, density);
    printf("sdl %g, density %g\n",sdl, density);
  }

}
