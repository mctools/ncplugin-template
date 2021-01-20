#include "NCSansHelper.hh"

bool NCPluginNamespace::calSDL(const NC::Info& info, double &scatLenDensity, double &numberDensity)
{
  scatLenDensity= 0.;
  numberDensity =0.;
  //calculate scattering length density from the dynamic info
  if(info.hasDynamicInfo()&&info.hasNumberDensity())
  {
    numberDensity = info.getNumberDensity();   // in atoms/Aa^3
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
    double perVolume = 1./strInfo.volume;//Aa^3

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

NC::Info::CustomSectionData::const_iterator NCPluginNamespace::findCustomLineIter(const NC::Info::CustomSectionData& data, const std::string& keyword, bool check)
{
  NC::Info::CustomSectionData::const_iterator it(data.end());
  for(auto line=data.begin();line!=data.end();++line)
  {
    if(line->at(0)==keyword)
    {
      it=line;
      break;
    }
  }
  if(check && it==data.end())
    NCRYSTAL_THROW2(BadInput,"findCustomLineIter can not find parameter " << keyword);
  return it;
}
