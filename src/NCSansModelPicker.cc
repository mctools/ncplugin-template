#include "NCSansModelPicker.hh"
#include  "NCSansHelper.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include <vector>

bool NCP::SansModelPicker::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::SansIsotropic NCP::SansModelPicker::createFromInfo( const NC::Info& info)
{
  auto iq = SansModelPicker(info);
  return SansIsotropic(iq.getQ(), iq.getI());
}

NCP::SansModelPicker::IqCalType NCP::SansModelPicker::getIqCalType(const NC::Info::CustomSectionData& data) const
{
  for(auto line:data)
  {
    if(line.at(0)=="DirectLoad")
      return IqCalType::kDirectLoad;
    else if(line.at(0)=="HardSphere")
      return IqCalType::kHardSphere;
  }
  return IqCalType::kUndefined;
}

void NCP::SansModelPicker::IqHardSphere(const NC::Info::CustomSectionData& data, const NC::Info& info)
{
  //radius
  auto it_r=findCustomLineIter(data, "radius");
  if(it_r->size()!=2)
    NCRYSTAL_THROW2( BadInput,"radius in the @CUSTOM_"<<pluginNameUpperCase()
                   <<" radius field should prove one parameter" );
  double radius(0.);
  if ( !NC::safe_str2dbl( it_r->at(1), radius ) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section radius" );

  //solvent
  NC::Info::CustomSectionData::const_iterator it_solvent=data.end();

  try {
    it_solvent=findCustomLineIter(data, "solvent");
  } catch (NC::Error::BadInput&e) {
    //I don't prevent no solvent
  }

  if(it_solvent!=data.end())
  {
    if(it_solvent->size()!=3)
      NCRYSTAL_THROW2( BadInput,"radius in the @CUSTOM_"<<pluginNameUpperCase()
                     <<" solvent field should prove two parameters" );
    m_solvantCfg = it_solvent->at(1);

    double solventfac (0.);
    if ( !NC::safe_str2dbl( it_solvent->at(2), solventfac ) )
      NCRYSTAL_THROW2( BadInput,"Invalid volume fraction specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" solvent " );
    m_solvantCfg += ";packfact="+std::to_string(solventfac);
    printf("m_solvantCfg %s\n", m_solvantCfg.c_str());
  }

  double sld(0), numden(0);
  calSDL(info, sld, numden);
  printf("sdl %g, number density %g\n",sld, numden);

  double R=radius;
  double R3=R*R*R;
  double V = 4./3.*NC::kPi*R3;
  double atomNumInSphere = V*numden;


  m_Q=NC::logspace(-6,10,1000);
  m_I.reserve(m_Q.size());

  for(double q:m_Q)
  {
    if(q<1e-5) // approximate by the limit at zero, fixme: this should be found automatically
      m_I.push_back(1.77777777777777777777777778*NC::kPi*NC::kPi*pow(radius,6)*sld*sld/atomNumInSphere);
    else
    {
      double P = 3* (sin(q*R) - q*R*cos(q*R))/(R3* q*q*q);
      m_I.push_back( V*V* sld* sld* P*P/atomNumInSphere);
    }
  }

}

void NCP::SansModelPicker::IqDirectLoad(const NC::Info::CustomSectionData& data)
{
  //Verify we have 3 lines and 2 vectors has identical number of elements
  // if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
  //   NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
  //                   <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  auto dataQ=findCustomLineIter(data, "Q");
  auto dataI=findCustomLineIter(data, "I");

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

NCP::SansModelPicker::SansModelPicker( const NC::Info& info )
:m_packfact(1.), m_volfact(1.)
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.


  // unsigned numSec =  info.countCustomSections( pluginNameUpperCase() ); fixme

    NC::Info::CustomSectionData data = info.getCustomSection( pluginNameUpperCase(), 0 );

    switch(getIqCalType(data)) {
      case kDirectLoad:
        IqDirectLoad(data);
        break;
      case kHardSphere:
        IqHardSphere(data,info);
        break;
      default :
        NCRYSTAL_THROW2(BadInput," @CUSTOM_"<<pluginNameUpperCase()<< " with undefined load method");
  }

}
